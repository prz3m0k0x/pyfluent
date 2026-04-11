from dataclasses import dataclass, field
from scipy.optimize import Bounds, LinearConstraint
import numpy as np
import pathlib, shutil, os, subprocess, json, csv, sys, traceback
from os import PathLike
import ansys.geometry.core as geometry
import ansys.fluent.core as pyfluent
from ansys.fluent.core.solver import (
    Species,
    Viscous,
    Mesh,
    Energy,
    Initialization,
    RunCalculation,
    ReportDefinitions,
)

@dataclass
class PSOConfig:
    
    POP_SIZE    : int = 15
    MAX_ITER    : int = 30
    N_PARAMS    : int = 4
    N_RESPONSES : int = 3

    w_init      : float = 0.8
    c1_init     : float = 0.6
    c2_init     : float = 0.6
    w_finish    : float = 0.5
    c1_finish   : float = 0.9
    c2_finish   : float = 0.9

    x_lb : list = field(default_factory=lambda: [500, 500, 500, 0.05])
    x_ub : list = field(default_factory=lambda: [1000, 2000, 2000, 0.15])

    constr_matrix : np.typing.ArrayLike = field(default_factory =  lambda:np.array([
                                                [1,1,1,0],
                                                [0,0,0,0],
                                                [0,0,0,0],
                                                [0,0,0,0]]))
    
    contrs_lb : np.typing.ArrayLike = field(default_factory =  lambda:np.array(np.array([
                                                -np.inf,
                                                -np.inf,
                                                -np.inf,
                                                -np.inf
                                                ])))

    constr_ub : np.typing.ArrayLike = field(default_factory =  lambda:np.array(
                                                np.array([4000,
                                                np.inf,
                                                np.inf,
                                                np.inf
                                                ])))

    bound_weights       : list = field(default_factory=lambda: [1, 1, 1, 1])
    objective_weights   : list= field(default_factory=lambda: [0.5, 0.25, 0.25])
    penalty_coeff       : float = 1.0

    @property
    def bounds(self) -> Bounds :
        return Bounds(self.x_lb, self.x_ub)
    
    @property
    def constraints(self) -> LinearConstraint :
        return LinearConstraint(self.constr_matrix, self.contrs_lb, self.constr_ub)
    
    @property
    def algorithm_velocity_parameters(self) -> np.typing.ArrayLike :
        return np.linspace(start = [self.w_init, self.c1_init, self.c2_init],
                                                         stop  = [self.w_finish, self.c1_finish, self.c2_finish],
                                                         num= self.MAX_ITER)
    
class SensitivityData:
    """
    Reads a sensitivity analysis .csv file, ranks design points, and
    provides the best ones as the initial swarm for PSO.
    """

    def __init__(self, filepath: PathLike, parameters: list, responses: list,
                 response_weights: np.ndarray, goals: list):
        self.filepath         = filepath
        self.parameters       = parameters
        self.responses        = responses
        self.response_weights = np.array(response_weights)
        self.goals            = goals

        self.par_dict, self.response_dict, self.headers, self.raw_data = \
            self._read_csv()

    def _read_csv(self):
        headers, data = [], []

        try:
            with open(self.filepath, mode="r") as f:
                reader = csv.reader(f, delimiter=";")
                for row_count, row in enumerate(reader):
                    if row_count == 0:
                        headers = row
                        continue
                    row_values = []
                    for col_idx, col_value in enumerate(row):
                        try:
                            if col_value.strip().lower() == "false":
                                print(f"Unfeasible solution - Row {row_count}!")
                                break
                            row_values.append(col_value)
                        except ValueError:
                            print(f"Error - Column {col_idx}, Row {row_count}!")
                    data.append(row_values)

            print(f"Read .csv with headers:\n{headers}")

        except Exception:
            raise ValueError(f"Cannot read file {os.path.basename(self.filepath)}!")

        data_T  = list(map(list, zip(*data)))
        id_list = list(range(len(data)))

        par_idx = [headers.index(p) for p in self.parameters]
        obj_idx = [headers.index(r) for r in self.responses]

        par_dict = []
        for idx in par_idx:
            col = [float(v) for v in data_T[idx]]
            par_dict.append({"header": headers[idx], "id": id_list, "data": col})

        response_dict = []
        for idx in obj_idx:
            col = [float(v) for v in data_T[idx]]
            response_dict.append({"header": headers[idx], "id": id_list, "data": col})

        return par_dict, response_dict, headers, data

    @staticmethod
    def _normalize(data: np.ndarray, weights: np.ndarray, goals: list) -> np.ndarray:
        """Returns NSOFV shape (n_points,)."""
        signs = np.array([-1 if g == "max" else 1 for g in goals])[:, np.newaxis]
        w     = weights[:, np.newaxis]
        means = np.mean(data, axis=1)[:, np.newaxis]
        return np.sum(w * signs * (data / means), axis=0)

    @staticmethod
    def _can_be_broadcast(*args) -> bool:
        try:
            np.broadcast(*args)
            return True
        except ValueError:
            return False

    def best_designs(self, n_designs: int):
        """
        Returns the n best design points ranked by weighted NSOFV.

        Returns
        -------
        best_params  : np.ndarray (N_PARAMS, n_designs)
        best_nsofv   : np.ndarray (n_designs,)
        best_ids     : np.ndarray of int (n_designs,)
        """
        data_matrix = np.array([r["data"] for r in self.response_dict])
        ids         = np.array(self.response_dict[0]["id"])

        nsofv    = self._normalize(data_matrix, self.response_weights, self.goals)

        id_nsofv = np.vstack((ids, nsofv))
        sort_idx = np.argsort(id_nsofv[1, :])
        sorted_  = id_nsofv[:, sort_idx]

        best_ids   = sorted_[0, :n_designs].astype(int)
        best_nsofv = sorted_[1, :n_designs]

        best_params = np.array([
            [param["data"][pos] for pos in best_ids]
            for param in self.par_dict
        ])                                           

        return best_params, best_nsofv, best_ids

    @property
    def response_means(self) -> np.ndarray:
        """Shape (N_RESPONSES, 1) — used as f_mean_init in PSOOptimizer."""
        data_matrix = np.array([r["data"] for r in self.response_dict])
        return np.mean(data_matrix, axis=1, keepdims=True)

class OptiRun:

    def __init__(self, cwd: pathlib.Path, run_path : pathlib.Path, solver_path: pathlib.Path,
                mixture_path: pathlib.Path, geometry_script: pathlib.Path) -> None:
        self.cwd             = cwd
        self.run_path        = pathlib.Path(run_path)
        self.solver_path     = pathlib.Path(solver_path)
        self.mixture_path    = pathlib.Path(mixture_path)
        self.geometry_script = pathlib.Path(geometry_script)
        self.opti_dir        = cwd / "opti_dir"
        self.py_exe          = sys.executable

        for label, path in [("solver",   self.solver_path),
                            ("mixture",  self.mixture_path),
                            ("geometry", self.geometry_script)]:
            if not path.exists():
                raise FileNotFoundError(f"{label} file not found: {path}")

        self.make_root_dir()

    def make_root_dir(self) -> None:
        parent_scripts = self.cwd / "parent-scripts"
        parent_scripts.mkdir(exist_ok=True)
        self.opti_dir.mkdir(exist_ok=True)

        for src in [self.solver_path, self.mixture_path, self.geometry_script, self.run_path]:
            shutil.copy(src, parent_scripts / src.name)

    @staticmethod
    def _conversion(so3: float, so2: float) -> float:
        M_so3, M_so2 = 80.06, 64.07
        return (so3 / M_so3) / ((so3 / M_so3) + (so2 / M_so2))

    @staticmethod
    def _extract_last_value(filepath: pathlib.Path) -> float:
        if filepath.exists():
            lines = [l for l in filepath.read_text().splitlines() if l.strip()]
            if lines:
                return float(lines[-1].split()[1])
        print(f"Warning: Could not read {filepath}")
        return 0.0

    @staticmethod
    def _check_convergence(filepath: pathlib.Path, max_iter: int) -> bool:

        if filepath.exists():
            lines = [l for l in filepath.read_text().splitlines() if l.strip()]
            if lines:
                converged = int(lines[-1].split()[0]) < max_iter
                status    = "converged" if converged else "did NOT converge"
                print(f"Case in {filepath.parent} {status}!")
                return converged
        print(f"Warning: Could not read {filepath}")
        return False

    def run(self, LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2,
            particle_dir: pathlib.Path) -> np.ndarray:
        particle_dir.mkdir(exist_ok=True)

        subprocess.run(
            [sys.executable, str(self.cwd / "run.py"),
            str(LengthCat1), str(LengthCool), str(LengthCat2), str(inlet_Y_SO2),
            str(particle_dir), str(self.geometry_script),
            str(self.mixture_path), str(self.solver_path)],
            check=True,
            cwd=particle_dir,
        )

        with open(particle_dir / "response.json") as f:
            result = json.load(f)

        return np.array([result["conversion"], result["so3"], result["so2"]])

    def mock_run(self, *args, **kwargs) -> np.ndarray:
        return np.array([np.random.uniform(0.7,  0.99),
                         np.random.uniform(0.05, 0.25),
                         np.random.uniform(0.01, 0.08)])

    def evaluate_swarm(self, particles: np.ndarray, iteration: int,
                       test_mode: bool = True) -> np.ndarray:
        n_responses, pop_size = 3, particles.shape[1]
        responses = np.zeros((n_responses, pop_size))
        call      = self.mock_run if test_mode else self.run

        for i in range(pop_size):
            particle_dir    = self.opti_dir / f"Iter{iteration:04d}_Particle{i:03d}"
            responses[:, i] = call(
                LengthCat1   = particles[0, i],
                LengthCool   = particles[1, i],
                LengthCat2   = particles[2, i],
                inlet_Y_SO2  = particles[3, i],
                particle_dir = particle_dir,
            )
        return responses
    
class Swarm:
    """Owns particle positions, velocities, personal and global bests."""

    def __init__(self, config: PSOConfig, initial_particles: np.ndarray,
                 initial_gains: np.ndarray, mode: str = "min"):
        self.config   = config
        self.mode     = mode
        self.particles = initial_particles.copy()
        self.velocity  = np.zeros_like(self.particles)

        self.pbest_positions = initial_particles.copy()
        self.pbest_gains     = initial_gains.copy()

        best_idx = np.argmin(initial_gains) if mode == "min" else np.argmax(initial_gains)
        self.gbest_position = initial_particles[:, best_idx].copy()
        self.gbest_gain     = initial_gains[best_idx]

    def step(self, w: float, c1: float, c2: float) -> None:
        """Updates velocity and positions. Violations handled by penalty only — no clamping."""
        n_params, pop_size = self.particles.shape
        r1 = np.random.uniform(0, 1, (n_params, pop_size))
        r2 = np.random.uniform(0, 1, (n_params, pop_size))

        self.velocity = (
            w  * self.velocity
          + c1 * r1 * (self.pbest_positions - self.particles)
          + c2 * r2 * (self.gbest_position[:, np.newaxis] - self.particles)
        )
        self.particles += self.velocity

    def update_bests(self, new_gains: np.ndarray) -> None:
        if self.mode == "min":
            improved = new_gains < self.pbest_gains
            is_new_best = new_gains[np.argmin(new_gains)] < self.gbest_gain
            best_idx    = np.argmin(new_gains)
        else:
            improved = new_gains > self.pbest_gains
            is_new_best = new_gains[np.argmax(new_gains)] > self.gbest_gain
            best_idx    = np.argmax(new_gains)

        self.pbest_positions[:, improved] = self.particles[:, improved]
        self.pbest_gains[improved]        = new_gains[improved]

        if is_new_best:
            self.gbest_position = self.particles[:, best_idx].copy()
            self.gbest_gain     = new_gains[best_idx]

class HistoryLogger:

    def __init__(self, config: PSOConfig):
        self.config            = config
        self.best_gain_history = []
        self.particle_history  = np.zeros((config.N_PARAMS, config.POP_SIZE, config.MAX_ITER))

    def log(self, iteration: int, particles: np.ndarray, best_gain: float) -> None:
        self.particle_history[:, :, iteration] = particles
        self.best_gain_history.append(best_gain)
        print(f"Iter {iteration+1:02d}  best_gain={best_gain:.6f}")

    def save(self, path: pathlib.Path) -> None:
        np.save(path / "particle_history.npy", self.particle_history)
        np.savetxt(path / "gain_history.csv",
                   np.array(self.best_gain_history), delimiter=";")

class SimulationSetup:

    # Mesh parameters
    GLOBAL_MIN: float = 60.
    GLOBAL_MAX: float = 240.
    LOCAL_SIZE: float = 120.

    # Simulation parameters:
    INLET_NAME          = "inlet"
    OUTLET_NAME         = "outlet"
    WALLS_NAME          = "walls"
    MIXTURE_MODEL       = "reaction-mixture"
    ZONE_FLUID1         = "fluid_zone1"
    ZONE_CAT1           = "fluid_catalyst1"
    ZONE_COOLING        = "fluid_cooling"
    ZONE_CAT2           = "fluid_catalyst2"
    ZONE_FLUID2         = "fluid_zone2"

    @property
    def ZoneAssignment(self):
        COOLING_ZONE        = field(default_factory=lambda: [self.ZONE_COOLING])
        CATALYST_ZONES      = field(default_factory=lambda: [self.ZONE_CAT1, self.ZONE_CAT2])

    VISCOUS_RESISTANCE  = 847936.476    #1/m^2
    INERTIAL_RESISTANCE = 283.104780 #1/m
    CATALYST_POROSITY = 0.6

    HEAT_EXCHANGER_VISCOUS_RESISTANCE = 1e4
    HEAT_EXCHANGER_INERTIAL_RESISTANCE = 100


    COOLING_HEAT_SINK   = -16000   #W/m^3

    # Inlet conditions
    INLET_VELOCITY      = 0.5
    INLET_TEMPERATURE   = 773.15 # K
    # INLET_Y_SO2         = inlet_Y_SO2       #SO2 mass fraction at inlet
    INLET_Y_O2          = 0.21 * (1 - INLET_Y_SO2)   #O2 in remaining air

    # Solver
    ITER_COUNT          = 500
    
    def __init__(self, cwd: pathlib.Path,
                 script_file: pathlib.Path,
                 mixture_file: pathlib.Path) -> None:
        self.cwd          = pathlib.Path(cwd).resolve()
        self.script_file  = pathlib.Path(script_file).resolve()
        self.mixture_file = pathlib.Path(mixture_file).resolve()

    def _run_geometry(self, script_args: dict) -> pathlib.Path:
        modeler = geometry.launch_modeler_with_spaceclaim()
        try:
            modeler.run_discovery_script_file(
                file_path   = str(self.script_file),
                script_args = script_args,
            )
            geometry_file = pathlib.Path(script_args["dir"]) / "geometry.pmdb"
            if not geometry_file.exists():
                raise FileNotFoundError(f"Geometry file not created: {geometry_file}")
            print(f"  [1/3] Geometry created: {geometry_file}")
            return geometry_file
        finally:
            modeler.exit()

    def _run_meshing(self, geometry_file: pathlib.Path,
                            particle_dir: pathlib.Path) -> pathlib.Path:

        mesh_path = particle_dir / "geometry.msh.h5"   # ← per-particle, not self.cwd


        meshing_session = pyfluent.launch_fluent(
            mode            = "meshing",
            precision       = "double",
            processor_count = 6,
            ui_mode         ="no_gui",
            cwd             = str(particle_dir),
            cleanup_on_exit = True
        )

        try:
            workflow = meshing_session.workflow  
            workflow.InitializeWorkflow(WorkflowType = "Watertight Geometry")
            
            workflow.TaskObject["Import Geometry"].Arguments = dict(FileName = str(geometry_file))
            workflow.TaskObject["Import Geometry"].Execute()

            local_sizing = meshing_session.workflow.TaskObject["Add Local Sizing"]
            local_sizing.Arguments.set_state(
                {
                    "AddChild"          : "yes",
                    "BOIControlName"    : "body_size1",
                    "BOIFaceLabelList"  : ["fluid_catalyst2", "fluid_catalyst1"],
                    "BOIExecution"      : "Body Size",
                    "BOISize"           : self.LOCAL_SIZE,
                }
            )
            local_sizing.AddChildAndUpdate()

            surface_mesh_gen = meshing_session.workflow.TaskObject["Generate the Surface Mesh"]
            surface_mesh_gen.Arguments.set_state(
                {"CFDSurfaceMeshControls": {"MaxSize": self.GLOBAL_MAX, "MinSize": self.GLOBAL_MIN}}
            )

            surface_mesh_gen.Execute()        


            workflow.TaskObject["Describe Geometry"].Arguments.set_state({
                "SetupType"     : "The geometry consists of only fluid regions with no voids",
                "WallToInternal": "Yes",
            })
            workflow.TaskObject["Describe Geometry"].Execute()


            workflow.TaskObject["Update Boundaries"].Execute()

            workflow.TaskObject["Update Regions"].Execute()

            workflow.TaskObject["Add Boundary Layers"].InsertCompoundChildTask()
            workflow.TaskObject["smooth-transition_1"].Arguments.set_state(
                {"NumberOfLayers": 10}
            )
            workflow.TaskObject["smooth-transition_1"].Execute()

            volume_mesh_gen = meshing_session.workflow.TaskObject["Generate the Volume Mesh"]
            volume_mesh_gen.Arguments.set_state(
                {
                    "VolumeFill"            : "poly-hexcore",
                    "VolumeMeshPreferences" : {
                        "CheckSelfProximity": "yes",
                        "ShowVolumeMeshPreferences": True,
                    },
                }
            )

            volume_mesh_gen.Execute()

            meshing_session.tui.file.write_mesh(mesh_path)

            meshing_session.tui.file.write_mesh(str(mesh_path))

        except Exception as e:
            print(f"Meshing failed for {geometry_file}")
            print("Type:", type(e).__name__)
            traceback.print_exc()
            raise

        finally:
            meshing_session.exit()

        print(f"  [2/3] Mesh created: {mesh_path}")
        return mesh_path
        
    def _run_solution(self, file_path: pathlib.Path,
                    particle_dir: pathlib.Path,       # ← renamed from cwd, used consistently
                    inlet_Y_SO2: float) -> None:

        INLET_Y_O2 = 0.21 * (1 - inlet_Y_SO2)          # ← computed locally, depends on particle

        so3_file = str(particle_dir / "so3-outlet-fraction.out")
        so2_file = str(particle_dir / "so2-outlet-fraction.out")
        case_path = str(particle_dir / "so2_reactor.cas.h5")

        solver_session = pyfluent.launch_fluent(
            mode            = pyfluent.FluentMode.SOLVER,
            cwd             = str(particle_dir),         # ← particle_dir, not self.cwd
            dimension       = 3,
            precision       = pyfluent.Precision.DOUBLE,
            processor_count = 8,
            cleanup_on_exit = True,
            ui_mode         = "no_gui",
        )

        try:
            solver_session.settings.file.read_mesh(file_name=str(file_path))  # ← str()

            mesh_obj = Mesh(solver_session, new_instance_name="mesh")
            mesh_obj.surfaces_list = mesh_obj.surfaces_list.allowed_values()
            mesh_obj.display()

            general = solver_session.settings.setup.general
            general.solver.time.set_state("steady")
            general.solver.type.set_state("pressure-based")

            viscous = Viscous(solver_session)
            viscous.model           = "k-epsilon"
            viscous.k_epsilon_model = "realizable"

            Energy(solver_session).enabled = True

            species = Species(solver_session)
            species.model.option = "species-transport"
            species.reactions.enable_volumetric_reactions = True

            try:
                solver_session.tui.define.materials.data_base.database_type(
                    "user-defined", str(self.mixture_file))
                solver_session.tui.define.materials.copy("mixture", self.MIXTURE_MODEL)
                print(f"  Loaded mixture from {self.mixture_file}")
            except Exception as e:
                print(f"  Warning: Failed to load mixture SCM: {e}")

            species.model.material = self.MIXTURE_MODEL

            fluid_zone = solver_session.settings.setup.cell_zone_conditions.fluid
            for fluid in fluid_zone.get_object_names():
                fluid_zone[fluid].reaction.react = False

            for zone in self.CATALYST_ZONES:
                fz = fluid_zone[zone]
                fz.porous_zone.porous = True
                fz.porous_zone.viscous_resistance.set_state(
                    [{"option": "value", "value": self.VISCOUS_RESISTANCE}] * 3)
                fz.porous_zone.inertial_resistance.set_state(
                    [{"option": "value", "value": self.INERTIAL_RESISTANCE}] * 3)
                fz.porous_zone.porosity.set_state(
                    {"option": "value", "value": self.CATALYST_POROSITY})
                fz.reaction.react = True

            fz_cool = fluid_zone[self.ZONE_COOLING]
            fz_cool.sources.set_state({"enable": True})
            fz_cool.sources.terms["energy"].resize(1)
            fz_cool.sources.terms.set_state({"energy": [{"option": "value"}]})
            fz_cool.sources.terms["energy"][0].set_state(self.COOLING_HEAT_SINK)
            fz_cool.porous_zone.porous = True
            fz_cool.porous_zone.viscous_resistance.set_state(
                [{"option": "value", "value": self.HEAT_EXCHANGER_VISCOUS_RESISTANCE}] * 3)
            fz_cool.porous_zone.inertial_resistance.set_state(   # ← was ZONE_COOLING bare name
                [{"option": "value", "value": self.HEAT_EXCHANGER_INERTIAL_RESISTANCE}] * 3)

            bcs = solver_session.settings.setup.boundary_conditions
            bcs.settings.physical_velocity_porous_formulation = True
            bcs.set_zone_type(zone_list=[self.INLET_NAME], new_type="velocity-inlet")
            inlet = bcs.velocity_inlet[self.INLET_NAME]
            inlet.momentum.velocity_magnitude                = self.INLET_VELOCITY
            inlet.thermal.temperature.value                  = self.INLET_TEMPERATURE
            inlet.species.species_mass_fraction["so2"].value = inlet_Y_SO2
            inlet.species.species_mass_fraction["o2"].value  = INLET_Y_O2  # ← local variable

            outlet = bcs.pressure_outlet[self.OUTLET_NAME]
            outlet.thermal.backflow_total_temperature.value = 700

            rep_defs = ReportDefinitions(solver_session)

            rep_defs.surface.create(name="so3-outlet-fraction")
            so3_rep = rep_defs.surface["so3-outlet-fraction"]
            so3_rep.report_type      = "surface-areaavg"
            so3_rep.field            = "so3"
            so3_rep.surface_names    = [self.OUTLET_NAME]
            so3_rep.output_parameter = True

            rep_defs.surface.create(name="so2-outlet-fraction")
            so2_rep = rep_defs.surface["so2-outlet-fraction"]
            so2_rep.report_type   = "surface-areaavg"
            so2_rep.field         = "so2"
            so2_rep.surface_names = [self.OUTLET_NAME]

            solver_session.settings.solution.monitor.report_files.create(name="so3-monitor")
            solver_session.settings.solution.monitor.report_files["so3-monitor"] = {
                "file_name": so3_file, "report_defs": ["so3-outlet-fraction"]}

            solver_session.settings.solution.monitor.report_files.create(name="so2-monitor")
            solver_session.settings.solution.monitor.report_files["so2-monitor"] = {
                "file_name": so2_file, "report_defs": ["so2-outlet-fraction"]}

            initialization = Initialization(solver_session)
            initialization.initialization_type = "hybrid"
            initialization.initialize()

            run_calc = RunCalculation(solver_session)
            run_calc.iter_count = self.ITER_COUNT
            run_calc.calculate()

            solver_session.settings.file.write_case_data(file_name=case_path)
            print(f"  [3/3] Solver complete. Case saved to {case_path}")

            # ← no return value — caller reads .out files directly

        except Exception as e:
            print(f"  Solver failed: {type(e).__name__}: {e}")
            traceback.print_exc()
            raise
        finally:
            solver_session.exit()          # ← always closes, regardless of success/failure
        
class PSOOptimizer:

    def __init__(self, config: PSOConfig, runner: OptiRun,
                 initial_particles: np.ndarray, f_mean_init: np.ndarray,
                 mode: str = "min"):
        self.config      = config
        self.runner      = runner
        self.f_mean_init = f_mean_init
        self.mode        = mode

        initial_responses = runner.evaluate_swarm(initial_particles, iteration=0)
        initial_gains     = self._compute_gains(initial_responses, initial_particles)

        self.swarm  = Swarm(config, initial_particles, initial_gains, mode=mode)
        self.logger = HistoryLogger(config)

    @classmethod
    def from_sensitivity(cls, config: PSOConfig, runner: OptiRun,
                         sens: SensitivityData, n_designs: int,
                         mode: str = "min") -> "PSOOptimizer":
        """Seeds the swarm from the best n_designs sensitivity analysis points."""
        best_params, _, _ = sens.best_designs(n_designs)
        f_mean_init       = sens.response_means
        return cls(config, runner, best_params, f_mean_init, mode=mode)

    @classmethod
    def from_random(cls, config: PSOConfig, runner: OptiRun,
                    mode: str = "min") -> "PSOOptimizer":
        """Seeds the swarm with uniformly random particles within bounds."""
        lb = np.array(config.x_lb)[:, np.newaxis]
        ub = np.array(config.x_ub)[:, np.newaxis]
        particles   = np.random.uniform(0, 1, (config.N_PARAMS, config.POP_SIZE)) * (ub - lb) + lb
        responses   = runner.evaluate_swarm(particles, iteration=0)
        f_mean_init = np.mean(responses, axis=1, keepdims=True)
        return cls(config, runner, particles, f_mean_init, mode=mode)


    def _compute_gains(self, responses: np.ndarray, particles: np.ndarray) -> np.ndarray:
        w   = np.array(self.config.objective_weights)[:, np.newaxis]
        obj = np.sum(w * (responses / self.f_mean_init), axis=0)
        return obj + self._penalty(particles)

    def _penalty(self, particles: np.ndarray) -> np.ndarray:
        bounds     = self.config.bounds
        constraint = self.config.constraints
        lb  = bounds.lb[:, np.newaxis]
        ub  = bounds.ub[:, np.newaxis]
        bw  = np.array(self.config.bound_weights)[:, np.newaxis]
        Ax  = constraint.A @ particles

        bound_viol  = np.sum(((lb > particles) | (particles > ub)).astype(int) * bw, axis=0)
        constr_viol = np.sum(
            (Ax > self.config.constr_ub[:, np.newaxis]).astype(int), axis=0)

        return self.config.penalty_coeff * (bound_viol + constr_viol)

    #main loop
    def run(self, test_mode: bool = True) -> Swarm:
        algo = self.config.algorithm_velocity_parameters

        for j in range(self.config.MAX_ITER):
            w, c1, c2 = algo[j]
            self.swarm.step(w, c1, c2)

            responses = self.runner.evaluate_swarm(self.swarm.particles,
                                                   iteration=j + 1,
                                                   test_mode=test_mode)
            new_gains = self._compute_gains(responses, self.swarm.particles)
            self.swarm.update_bests(new_gains)
            self.logger.log(j, self.swarm.particles, self.swarm.gbest_gain)

        return self.swarm
    
if __name__ == "__main__":
    config = PSOConfig()
    runner = OptiRun(
        cwd             = pathlib.Path.cwd(),
        run_path        = pathlib.Path(r"C:\Users\hp\Desktop\Me\Studia\CFD\pyfluent\pyfluent\so2\run.py"),
        solver_path     = pathlib.Path(r"C:\Users\hp\Desktop\Me\Studia\CFD\pyfluent\pyfluent\so2\solution.py"),
        mixture_path    = pathlib.Path(r"C:\Users\hp\Desktop\Me\Studia\CFD\pyfluent\pyfluent\so2\reaction-mixture.scm"),
        geometry_script = pathlib.Path(r"C:\Users\hp\Desktop\Me\Studia\CFD\pyfluent\pyfluent\so2\geometry_script_3d_python_v.scscript"),
    )
    (pathlib.Path.cwd() / "opti_dir").mkdir(exist_ok=True)

    sens = SensitivityData(
        filepath         = r"C:\Users\hp\Desktop\Me\Studia\CFD\pyfluent\pyfluent\so2\designpoints.csv",
        parameters       = ["inlet_Y_SO2", "LengthCat2", "LengthCool", "LengthCat1"],
        responses        = ["conversion", "so2", "so3"],
        response_weights = np.array([0.4, 0.3, 0.3]),
        goals            = ["max", "min", "min"],
    )

    print("response_means shape:", sens.response_means.shape)
    best_params, best_nsofv, best_ids = sens.best_designs(n_designs=15)
    print("best_params shape:", best_params.shape)
    print("best_ids:", best_ids)
    print("best_nsofv:", best_nsofv.round(4))

    optimizer = PSOOptimizer.from_sensitivity(config, runner, sens,
                                              n_designs=15, mode="min")
    optimizer.run(test_mode=False)