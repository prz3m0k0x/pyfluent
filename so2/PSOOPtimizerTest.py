from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import Bounds, LinearConstraint
from scipy.spatial import distance_matrix
import numpy as np
import numpy.typing as typing
import pathlib, shutil, os, json, csv, sys, traceback, time, math
from typing import Callable, List
from os import PathLike
import ansys.geometry.core as geometry
from ansys.geometry.core import launch_modeler_with_spaceclaim
from ansys.geometry.core.math import (
    Plane, Point2D, Point3D,
    UNITVECTOR3D_X, UNITVECTOR3D_Y, UNITVECTOR3D_Z,
)
from ansys.geometry.core.designer.design import DesignFileFormat
from ansys.geometry.core.misc import UNITS, Distance
from ansys.geometry.core.sketch import Sketch
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


def TestFunction(particles: np.ndarray):
    """
    Multi-objective test function.
    Objective 1 : quadratic form, minumim at x=(0,1), f=0
    Objective 2 : Rastrigin function — global min at x=(0,0), f=0

    particles : (n_params, N_PARTICLES)
    returns   : y1 (N_PARTICLES,), y2 (N_PARTICLES,)
    """
    x1 = particles[0, :]
    x2 = particles[1, :]

    #quadratic
    y1 = (x1-3)**2 + (x2-3)**2
    y2 = (x1+5)**2 + (x2+5)**2
    #  Rastrigin
    # A  = 10
    # y2 = (2 * A
    #        + x1**2 - A * np.cos(2 * np.pi * x1)
    #        + x2**2 - A * np.cos(2 * np.pi * x2))

    return y1, y1


@dataclass
class PSOConfig:
    h_factor    : int = 60
    max_iter    : int = 40
    n_params    : int = 2
    n_responses : int = 2
    t_neighbors : int = 10

    w_init      : float = 0.7
    c1_init     : float = 0.6
    c2_init     : float = 0.5
    w_finish    : float = 0.4
    c1_finish   : float = 0.8
    c2_finish   : float = 0.6
    v_max_factor: float = 0.2

    x_lb : list = field(default_factory=lambda: [])
    x_ub : list = field(default_factory=lambda: [])

    constr_matrix : list = field(default_factory=lambda: [])
    constr_lb : list = field(default_factory=lambda: [])
    constr_ub : list = field(default_factory=lambda: [])
    simulation_factory : Callable = field(default=None, repr=True)
    


    def __post_init__(self) -> None:
        self.x_lb = np.full(self.n_params, -np.inf, dtype=float) if len(self.x_lb) == 0 else np.array(self.x_lb, dtype=float)
        self.x_ub = np.full(self.n_params,  np.inf, dtype=float) if len(self.x_ub) == 0 else np.array(self.x_ub, dtype=float)

        self.constr_lb = np.array(self.constr_lb, dtype=float) if len(self.constr_lb) > 0 else np.array([], dtype=float)
        self.constr_ub = np.array(self.constr_ub, dtype=float) if len(self.constr_ub) > 0 else np.array([], dtype=float)
        self.constr_matrix = np.array(self.constr_matrix, dtype=float) if len(self.constr_matrix) > 0 else np.empty((0, self.n_params), dtype=float)

        if self.x_lb.shape != (self.n_params,):
            raise ValueError("x_lb must have length n_params")
        if self.x_ub.shape != (self.n_params,):
            raise ValueError("x_ub must have length n_params")

        if self.constr_matrix.shape[0] > 0:
            if self.constr_matrix.shape[1] != self.n_params:
                raise ValueError("constr_matrix must have n_params columns")
            n_constr = self.constr_matrix.shape[0]
            if self.constr_lb.size not in (0, n_constr):
                raise ValueError("constr_lb must have one value per constraint")
            if self.constr_ub.size not in (0, n_constr):
                raise ValueError("constr_ub must have one value per constraint")

    def POP_SIZE(self) -> int:
        if self.h_factor <= self.n_responses:
            raise ValueError("h_factor must be greater than number of responses!")
        return math.comb(self.h_factor + self.n_responses - 1, self.n_responses - 1)
    

    @staticmethod
    def das_dennis_weights(m: int, H: int) -> np.ndarray:
        """
        Generate weight matrix using Das Dennis method.
        m : number of objective functions
        H : partition parameter (divisions per axis)
        Returns: np.ndarray of shape (N, m), each row is a weight vector
        """
        def enumerate_weights(m, H, current=[]):
            if m == 1:
                yield current + [H]
            else:
                for i in range(H + 1):
                    yield from enumerate_weights(m - 1, H - i, current + [i])
        
        weights = np.array(list(enumerate_weights(m, H)), dtype=float) / H
        return weights
    
    @property
    def algorithm_velocity_parameters(self) -> np.typing.ArrayLike :
        return np.linspace(start = [self.w_init, self.c1_init, self.c2_init],
                                                         stop  = [self.w_finish, self.c1_finish, self.c2_finish],
                                                         num= self.max_iter)
@dataclass
class GeometrySetup:
    zones       : List[str]   = field(default_factory=lambda: ["In_zone", "LengthCat1", "LengthCool", "LengthCat2", "Out_zone"])
    inlet       : str = "inlet"
    outlet      : str = "outlet"
    walls       : str = "walls"
    lengths     : List[float] = field(default_factory=lambda: [500.0, 600.0, 1200.0, 1200.0, 2500.0])
    name        : str         = "geometry"
    radius      : float       = 1250.
    component   : str         = "reactor"
@dataclass
class MeshSetup:

    PROCESSOR_COUNT  : int   = 6
    PRECISION        : str   = "double"
    GUI              : str   = "no_gui"
    CLEANUP_ON_EXIT  : bool  = True

    WORKFLOWTYPE     : str   = "Watertight Geometry"

    GLOBAL_MIN       : float = 60.
    GLOBAL_MAX       : float = 240.

    N_LAYERS         : int   = 10
    INFLATION_METHOD : str   = "smooth-transition"
    INFLATION_SCOPE  : List[str] = field(default_factory=lambda: ["walls"])

    catalyst_zones   : List[str] = field(default_factory=lambda: ["reactor-lengthcat1", "reactor-lengthcat2"])
    cooling_zones    : List[str] = field(default_factory=lambda: ["reactor-lengthcool"])
    catalyst_size    : float = 120.
    cooling_size     : float = 180.

    sizing_configs   : List[dict] = field(init=False)
    boundary_layer_configs : List[dict] = field(init=False)
    global_sizing_configs : List[dict] = field(init=False)

    def __post_init__(self):
        self.sizing_configs = [
            self.SizingDictionary(
                boi_face_label_list = self.catalyst_zones,
                boiexecution        = "Body Size",
                local_size          = self.catalyst_size,
                boi_name            = "catalyst_size",
            ),
            self.SizingDictionary(
                boi_face_label_list = self.cooling_zones,
                boiexecution        = "Body Size",
                local_size          = self.cooling_size,
                boi_name            = "cooling_size",
            ),
        ]

        self.boundary_layer_configs = [
            self.BoundaryLayerDictionary(
                n_layers         = self.N_LAYERS,
                offset_method    = self.INFLATION_METHOD,
                bl_label_list    = self.INFLATION_SCOPE,
                bl_control_name  = "bl_walls",
            ),
        ]

        self.global_sizing_configs = self.GlobalSizingDictionary(
                                        volume_fill= "poly-hexcore",
                                        check_self_proximity= "yes")

    @staticmethod
    def BoundaryLayerDictionary(n_layers        : int,
                                offset_method   : str,
                                bl_label_list   : List[str]  = None,
                                zone_selection  : List[str]  = None,
                                region_scope    : List[str]  = None,
                                bl_control_name : str        = "bl_1",
                                transition_ratio: float      = 0.272,
                                addchild        : str        = "yes") -> dict:
        d = {
            "AddChild"         : addchild,
            "BLControlName"    : bl_control_name,
            "OffsetMethodType" : offset_method,
            "NumberOfLayers"   : n_layers,
            "TransitionRatio"  : transition_ratio,
        }
        if bl_label_list  is not None: d["BlLabelList"]       = bl_label_list
        if zone_selection is not None: d["ZoneSelectionList"] = zone_selection
        if region_scope   is not None: d["RegionScope"]       = region_scope
        return d
    
    @staticmethod
    def SizingDictionary(boi_face_label_list : List[str],
                         local_size          : float,
                         boiexecution        : str = "Body Size",
                         addchild            : str = "yes",
                         boi_name            : str = "body_size1") -> dict:
        return {
            "AddChild"         : addchild,
            "BOIControlName"   : boi_name,
            "BOIFaceLabelList" : boi_face_label_list,
            "BOIExecution"     : boiexecution,
            "BOISize"          : local_size,
        }
    
    @staticmethod
    def GlobalSizingDictionary(
        volume_fill : str = "poly-hexcore",
        check_self_proximity :str = "yes",):

        return {
            "VolumeFill" : volume_fill,
            "VolumeMeshPreferences" : {"CheckSelfProximity" : check_self_proximity,
                                       "ShowVolumeMeshPreferences" : True}
        }

@dataclass
class SimulationSetup:
    mixture_file        : str   = r"C:\Users\hp\...\reaction-mixture.scm"
    MIXTURE_MODEL       : str      = "reaction-mixture"
    PROCESSOR_COUNT     : int   = 6


    VISCOUS_RESISTANCE                  : float = 847936.476
    INERTIAL_RESISTANCE                 : float = 283.104780
    CATALYST_POROSITY                   : float = 0.6

    HEAT_EXCHANGER_VISCOUS_RESISTANCE   : float = 1e4
    HEAT_EXCHANGER_INERTIAL_RESISTANCE  : float = 100.
    COOLING_HEAT_SINK                   : float = -16000.

    INLET_VELOCITY              : float = 0.5
    INLET_TEMPERATURE           : float = 773.15
    OUTLET_BACKFLOW_TEMPERATURE : float = 1073.15
    ITER_COUNT                  : int   = 100

    RESPONSES           : List[str]     = field(default_factory= lambda: ["so3-outlet-fraction", "so2-outlet-fraction"])
    REPSONSES_FIELD     : List[str]     = field(default_factory= lambda: ["so3", "so2"]) #actual names of fields from fluent, like pressure, velocity-magnitude etc.
    RESPONSES_SCOPE     : List[str]     = field(default_factory= lambda: ["surface-areaavg", "surface-areaavg"])
    RESPONSES_SURFACE   : List[str]     = field(default_factory= lambda: ["outlet", "outlet"])
@dataclass

class Swarm:
    """Owns particle positions, velocities, personal and global bests."""
    def __init__(self, config : PSOConfig, initial_particles: np.ndarray,
                initial_gains: np.ndarray):
        # initial_particles : (n_params, N_PARTICLES)
        # initial_gains     : (n_responses, N_PARTICLES)

        self.config   = config
        self.particles = initial_particles.copy()
        self.velocity  = np.zeros_like(self.particles)
        self.v_max     = config.v_max_factor * (
                            np.array(config.x_ub) - np.array(config.x_lb)
                        )[:, np.newaxis] 

        self.weights = config.das_dennis_weights(
            m=initial_gains.shape[0], H=config.h_factor
        )                                                                   # (N_PARTICLES, n_responses)

        #Neighbourhood
        diff           = self.weights[:, None, :] - self.weights[None, :, :]  # (N, N, M)
        dist_matrix    = np.linalg.norm(diff, axis=-1)                        # (N, N)
        T              = min(config.t_neighbors, dist_matrix.shape[0])
        self.neighbors = np.argsort(dist_matrix, axis=1)[:, :T]               # (N_PARTICLES, T)


        self.z_ref = np.min(initial_gains, axis=1)                            # (n_responses,)
        self.z_nad = np.max(initial_gains, axis=1)

        self.pbest_positions = initial_particles.copy()                       # (n_params, N_PARTICLES)
        self.pbest_gains     = initial_gains.copy()                           # (n_responses, N_PARTICLES)
        self.pbest_scalars   = self._tcheby_scalars(initial_gains, self.particles)            # (N_PARTICLES,)

        self.gbest_positions = np.zeros_like(initial_particles)               # (n_params, N_PARTICLES)
        self.gbest_gains     = np.zeros_like(initial_gains)                   # (n_responses, N_PARTICLES)
        self._update_neighborhood_bests()

    @property
    def global_best_position(self) -> np.ndarray:
        """Single best particle across the whole swarm."""
        return self.pbest_positions[:, np.argmin(self.pbest_scalars)].copy()

    def _penalty(self, particles: np.ndarray) -> np.ndarray:
        lb  = np.array(self.config.x_lb)[:, np.newaxis]
        ub  = np.array(self.config.x_ub)[:, np.newaxis]

        #Quadratic bound violation penalty function handler
        lb_viol = np.maximum(0, lb - particles) ** 2          # (n_params, N_PARTICLES)
        ub_viol = np.maximum(0, particles - ub) ** 2
        bound_penalty = np.sum(lb_viol + ub_viol, axis=0)     # (N_PARTICLES,)

        Ax      = self.config.constr_matrix @ particles        # (N_CONSTR, N_PARTICLES)
        clb     = self.config.constr_lb[:, np.newaxis]
        cub     = self.config.constr_ub[:, np.newaxis]
        c_viol  = (np.maximum(0, clb - Ax) ** 2
                + np.maximum(0, Ax  - cub) ** 2)
        constr_penalty = np.sum(c_viol, axis=0)                # (N_PARTICLES,)

        return (bound_penalty + constr_penalty)
    
    def _update_neighborhood_bests(self) -> None:
        """
        For each particle i, find the neighbour with the lowest pbest_scalar
        and set it as that particle's gbest.
        """
        for i, neighbors in enumerate(self.neighbors):
            best_in_neighborhood   = neighbors[np.argmin(self.pbest_scalars[neighbors])]
            self.gbest_positions[:, i] = self.pbest_positions[:, best_in_neighborhood]
            self.gbest_gains[:, i]     = self.pbest_gains[:, best_in_neighborhood]

    def _tcheby_scalars(self, gains: np.ndarray,
                        particles: np.ndarray) -> np.ndarray:
        scale      = np.where(np.abs(self.z_nad - self.z_ref) > 1e-8,
                            self.z_nad - self.z_ref, 1.0)

        deviations = np.abs(gains - self.z_ref[:, np.newaxis]) / scale[:, np.newaxis]
        weighted   = self.weights.T * deviations
        tcheby     = np.max(weighted, axis=0)

        return tcheby + self._penalty(particles)


    def step(self, w: float, c1: float, c2: float) -> None:
        n_params, pop_size = self.particles.shape
        r1 = np.random.uniform(0, 1, (n_params, pop_size))
        r2 = np.random.uniform(0, 1, (n_params, pop_size))

        self.velocity = (
            w  * self.velocity
        + c1 * r1 * (self.pbest_positions - self.particles)
        + c2 * r2 * (self.gbest_positions - self.particles)
        )
        self.velocity  = np.clip(self.velocity, -self.v_max, self.v_max)
        self.particles += self.velocity

    def update_bests(self, new_gains: np.ndarray) -> None:

        self.z_ref = np.minimum(self.z_ref, np.min(new_gains, axis=1))
        self.z_nad = np.maximum(self.z_nad, np.max(new_gains, axis=1))


        new_scalars = self._tcheby_scalars(new_gains, self.particles)

        improved = new_scalars < self.pbest_scalars
        self.pbest_positions[:, improved] = self.particles[:, improved]
        self.pbest_gains[:, improved]     = new_gains[:, improved]
        self.pbest_scalars[improved]      = new_scalars[improved]


        self._update_neighborhood_bests()

class HistoryLogger:

    def __init__(self, config: PSOConfig):
        self.config            = config
        self.best_gain_history = []
        self.particle_history  = np.zeros((config.n_params, config.POP_SIZE(), config.max_iter))

    def log(self, iteration: int, particles: np.ndarray, best_gain: float) -> None:
        self.particle_history[:, :, iteration] = particles
        self.best_gain_history.append(best_gain)
        print(f"Iter {iteration+1:02d}  best_gain={best_gain:.6f}")

    def save(self, path: pathlib.Path) -> None:
        np.save(path / "particle_history.npy", self.particle_history)
        np.savetxt(path / "gain_history.csv",
                   np.array(self.best_gain_history), delimiter=";")
   
class PSOOptimizer:

    def __init__(self, config: PSOConfig, initial_particles: np.ndarray):
        self.config = config

        initial_gains = self._evaluate(initial_particles)               # (n_responses, N_PARTICLES)

        self.swarm  = Swarm(config, initial_particles, initial_gains)
        self.logger = HistoryLogger(config)

    @classmethod
    def from_random(cls, config: PSOConfig) -> "PSOOptimizer":
        """Seeds the swarm with uniformly random particles within bounds."""
        lb        = np.array(config.x_lb)[:, np.newaxis]
        ub        = np.array(config.x_ub)[:, np.newaxis]
        pop_size  = config.POP_SIZE()
        particles = np.random.uniform(0, 1, (config.n_params, pop_size)) * (ub - lb) + lb
        return cls(config, particles)


    def _evaluate(self, particles: np.ndarray, iteration: int = 0) -> np.ndarray:
        if self.config.simulation_factory is None:
            y1, y2 = TestFunction(particles)
            return np.vstack([y1, y2])

        n_particles = particles.shape[1]
        responses   = np.zeros((self.config.n_responses, n_particles))

        for i in range(n_particles):
            LengthCat1  = float(particles[0, i])
            LengthCool  = float(particles[1, i])

            so3, so2, conversion = self.config.simulation_factory(
                LengthCat1, LengthCool,
                iteration=iteration, particle=i,
            )
            responses[0, i] = -conversion
            responses[1, i] =  so2

        return responses

    def run(self) -> Swarm:
        algo             = self.config.algorithm_velocity_parameters
        particle_history = [self.swarm.particles.copy()]
        gbest_history    = [self.swarm.global_best_position.copy()]

        for j in range(self.config.max_iter):
            w, c1, c2 = algo[j]
            self.swarm.step(w, c1, c2)

            new_gains = self._evaluate(self.swarm.particles, iteration=j)
            self.swarm.update_bests(new_gains)

            self.logger.log(j, self.swarm.particles,
                            float(np.min(self.swarm.pbest_scalars)))
            particle_history.append(self.swarm.particles.copy())
            gbest_history.append(self.swarm.global_best_position.copy())

        visualize_pso(self.config, particle_history, gbest_history)
        return self.swarm
    
class VisualizePSO:
    def __init__(self, config : PSOConfig, history : HistoryLogger, swarm : Swarm, optimizer : PSOOptimizer):
        gbest_history = history.best_gain_history
        particle_history = history.particle_history
        self.optimizer = optimizer
        self.config = config
        self.swarm = swarm
            
    def visualize_pso(config : PSOConfig, particle_history, gbest_history,
                    trail_len: int = 15):
        """
        Animate PSO particle movement with a fading historical trail.

        trail_len : how many past iterations to show as fading ghost positions
        """
    def grid_gen(self, axes : typing.ArrayLike = (0,1), fixed_values : typing.ArrayLike = None, n_points : int = 200):

        if fixed_values == None:
            fixed_values = 0.5 * (np.array(self.config.x_lb) + np.array(self.config.x_ub))
        
        a,b = axes
        xa = np.linspace(self.config.x_lb[a], self.config.x_lb[a], n_points)
        xb = np.linspace(self.config.x_lb[b], self.config.x_lb[b], n_points)

        Xa, Xb = np.meshgrid(xa, xb)

        particles = np.tile(fixed_values[:, None], (1, Xa.size))
        particles[a, :] = Xa.ravel()
        particles[b, :] = Xb.ravel()
        
        return Xa, Xb, particles
    
    def tcheby_scalar_field(self, gains: np.ndarray, mode="min", weight_index=0):
        zref = self.swarm.zref
        znad = self.swarm.znad
        scale = np.where(np.abs(znad - zref) > 1e-12, znad - zref, 1.0)

        deviations = np.abs(gains - zref[:, None]) / scale[:, None]

        if mode == "single":
            w = self.swarm.weights[weight_index]
            return np.max(w[:, None] * deviations, axis=0)

        if mode == "min":
            scalars = [
                np.max(w[:, None] * deviations, axis=0)
                for w in self.swarm.weights
            ]
            return np.min(np.vstack(scalars), axis=0)

        raise ValueError("mode must be 'single' or 'min'")
    
    def evaluate_particles(self, particles : typing.ArrayLike):
        return self.optimizer._evaluate(particles)
    
    def plotter(self, plot_type : str, axes : typing.ArrayLike, f_vals : typing.ArrayLike = None, test_func : bool = True, n_points : int = 400):
        """Plot types supported:
        \"tscheby\" - creates tshebyshev scalar plot in chosen coordinate system (one parameter or two parameter) from all function outcomes
        \"func\" - creates 3d plot for chosen function in given coordinate system (one parameter or two parameter)
        \"pareto\" - creates 2d pareto front of best designs in given functional cooridnate system (f1, f2, ...)
        """
        if plot_type == "tscheby":
            grid = self.GridGen()
            values = TestFunction(grid)
            z_ref   = np.min(values, axis=1)
            Z = np.maximum(1 / float(self.config.n_responses * np.abs(values - z_ref))).reshape(grid[0].shape)
            fig, ax = plt.subplots(figsize=(9, 7))
            fig.patch.set_facecolor('#0d0d1a')
            ax.set_facecolor('#0d0d1a')

            X1 = grid[params[0]]
            X2 = grid[params[1]]

            hmap = ax.contourf(X1, X2, Z, levels=60, cmap='plasma', alpha=0.85)
            ax.contour(X1, X2, Z, levels=18, colors='white', alpha=0.10, linewidths=0.4)
            cbar = fig.colorbar(hmap, ax=ax, pad=0.02)
            cbar.set_label('Tchebycheff scalar', color='white', fontsize=10)
            cbar.ax.yaxis.set_tick_params(color='white')
            plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
        



class GeometryScript:
    def __init__(self, GeometrySetup : GeometrySetup, cwd : PathLike):
        self.geometry_setup = GeometrySetup
        self.cwd = pathlib.Path(cwd)
        self.filename = os.path.join(cwd, self.geometry_setup.name + ".pmdb")
        self.lengths = self.geometry_setup.lengths
        self.radius = Distance(self.geometry_setup.radius, unit= UNITS.mm)
        self.component = self.geometry_setup.component
        self.zones = self.geometry_setup.zones

    def _run_geometry(self):
        modeler = launch_modeler_with_spaceclaim()
        design  = modeler.create_design(self.geometry_setup.name)
        component = design.add_component(self.component)

        x_offset = 0.0
        bodies = []
        zones = self.geometry_setup.zones
        lengths = self.lengths
        radius = self.radius

        for zone, zone_length in zip(zones, lengths):
            plane = Plane(
                origin=Point3D([x_offset, 0, 0], unit=UNITS.mm),
                direction_x=UNITVECTOR3D_Y,
                direction_y=UNITVECTOR3D_Z,
            )
            sketch = Sketch(plane=plane)
            sketch.circle(
                center=Point2D([0, 0], unit=UNITS.mm),
                radius=radius,
            )
            body = component.extrude_sketch(
                name=zone,
                sketch=sketch,
                distance=Distance(zone_length, unit=UNITS.mm),
            )
            bodies.append(body)
            x_offset += zone_length

        inlet_face  = bodies[0].faces[0]
        outlet_face = bodies[-1].faces[1]
        wall_faces  = [face for body in bodies for j, face in enumerate(body.faces)
                    if j == 2]

        design.create_named_selection(self.geometry_setup.inlet,  faces=[inlet_face])
        design.create_named_selection(self.geometry_setup.outlet, faces=[outlet_face])
        design.create_named_selection(self.geometry_setup.walls,  faces=wall_faces)

        design.download(
            file_location = self.cwd / (self.geometry_setup.name + ".pmdb"),
            format        = DesignFileFormat.PMDB
        )
        modeler.close()
        return self.filename

class MeshScript:
    def __init__(self, GeometryScript : GeometryScript, MeshSetup : MeshSetup, geometry_file : PathLike) -> None:
        self.geometry_script = GeometryScript
        self.cwd = pathlib.Path(self.geometry_script.cwd)
        self.geometry_file = geometry_file
        self.setup = MeshSetup

    def _run_meshing(self) -> PathLike:
        cwd = pathlib.Path(self.cwd)
        mesh_path = pathlib.Path(cwd) / "geometry.msh.h5"


        meshing_session = pyfluent.launch_fluent(
            mode            = "meshing",
            precision       = self.setup.PRECISION,
            processor_count = self.setup.PROCESSOR_COUNT,
            ui_mode         = self.setup.GUI,
            cwd             = pathlib.Path(cwd),
            cleanup_on_exit = self.setup.CLEANUP_ON_EXIT
        )

        try:
            workflow = meshing_session.workflow
            workflow.InitializeWorkflow(WorkflowType = self.setup.WORKFLOWTYPE)
            
            workflow.TaskObject["Import Geometry"].Arguments = dict(FileName = str(self.geometry_file))
            workflow.TaskObject["Import Geometry"].Execute()

            
            for sizing in self.setup.sizing_configs:
                local_sizing = meshing_session.workflow.TaskObject["Add Local Sizing"]
                local_sizing.Arguments.set_state(sizing)
                local_sizing.AddChildAndUpdate()

            surface_mesh_gen = meshing_session.workflow.TaskObject["Generate the Surface Mesh"]
            surface_mesh_gen.Arguments.set_state(
                {"CFDSurfaceMeshControls": {"MaxSize": self.setup.GLOBAL_MAX, "MinSize": self.setup.GLOBAL_MIN}}
            )

            surface_mesh_gen.Execute()        


            workflow.TaskObject["Describe Geometry"].Arguments.set_state({
                "SetupType"     : "The geometry consists of only fluid regions with no voids",
                "WallToInternal": "Yes",
            })
            workflow.TaskObject["Describe Geometry"].Execute()


            workflow.TaskObject["Update Boundaries"].Execute()

            workflow.TaskObject["Update Regions"].Execute()

            for bl in self.setup.boundary_layer_configs:
                workflow.TaskObject["Add Boundary Layers"].InsertCompoundChildTask()

                child_name = f"{self.setup.INFLATION_METHOD}_{len(self.setup.boundary_layer_configs)}"
                workflow.TaskObject[child_name].Arguments.set_state(bl)
                workflow.TaskObject[child_name].Execute()


            volume_mesh_gen = meshing_session.workflow.TaskObject["Generate the Volume Mesh"]
            volume_mesh_gen.Arguments.set_state(
                self.setup.global_sizing_configs
                )

            volume_mesh_gen.Execute()

            meshing_session.tui.file.write_mesh(str(mesh_path))

        except Exception as e:
            print(f"Meshing failed for {self.geometry_file}")
            print("Type:", type(e).__name__)
            traceback.print_exc()
            raise

        finally:
            meshing_session.exit()

        print(f"  [2/3] Mesh created: {mesh_path}")

        return mesh_path

class SimulationScript:
    
    def __init__(self,
                 mesh_file : PathLike,
                 SimulationSetup: SimulationSetup,
                 MeshScript: MeshScript,
                 GeometrySetup : GeometrySetup) -> None:
        
        self.cwd          = pathlib.Path(MeshScript.cwd)
        self.setup        = SimulationSetup
        self.mesh_file    = mesh_file
        self.geometry_setup = GeometrySetup
        self.INLET_NAME = self.geometry_setup.inlet
        self.OUTLET_NAME = self.geometry_setup.outlet
        self.WALLS_NAME = self.geometry_setup.walls
        self.MIXTURE_MODEL = self.setup.MIXTURE_MODEL


        self.CATALYST_ZONES = [zone for zone in self.geometry_setup.zones if "cat" in str.lower(zone)]
        self.ZONE_COOLING = [zone for zone in self.geometry_setup.zones if "cool" in str.lower(zone)]

    def _run_simulation(self) -> tuple:
        


        cwd = pathlib.Path(self.cwd)
        inlet_y_o2 = 0.21 * (1 - 0.07)
        case_path = str(cwd / "so2_reactor.cas.h5")

        try:

            solver_session = pyfluent.launch_fluent(
                mode            = pyfluent.FluentMode.SOLVER,
                cwd             = pathlib.Path(cwd),      
                dimension       = 3,
                precision       = pyfluent.Precision.DOUBLE,
                processor_count = self.setup.PROCESSOR_COUNT,
                cleanup_on_exit = True,
                ui_mode         = "no_gui",                
            )

            solver_session.settings.file.read_mesh(file_name=str(self.mesh_file))  # ← str()

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
                    "user-defined", str(self.setup.mixture_file))
                solver_session.tui.define.materials.copy("mixture", self.setup.MIXTURE_MODEL)
                print(f"  Loaded mixture from {self.setup.mixture_file}")
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
                    [{"option": "value", "value": self.setup.VISCOUS_RESISTANCE}] * 3)
                fz.porous_zone.inertial_resistance.set_state(
                    [{"option": "value", "value": self.setup.INERTIAL_RESISTANCE}] * 3)
                fz.porous_zone.porosity.set_state(
                    {"option": "value", "value": self.setup.CATALYST_POROSITY})
                fz.reaction.react = True

            fz_cool = fluid_zone[self.ZONE_COOLING[0]]
            fz_cool.sources.set_state({"enable": True})
            fz_cool.sources.terms["energy"].resize(1)
            fz_cool.sources.terms.set_state({"energy": [{"option": "value"}]})
            fz_cool.sources.terms["energy"][0].set_state(self.setup.COOLING_HEAT_SINK)
            fz_cool.porous_zone.porous = True
            fz_cool.porous_zone.viscous_resistance.set_state(
                [{"option": "value", "value": self.setup.HEAT_EXCHANGER_VISCOUS_RESISTANCE}] * 3)
            fz_cool.porous_zone.inertial_resistance.set_state(
                [{"option": "value", "value": self.setup.HEAT_EXCHANGER_INERTIAL_RESISTANCE}] * 3)

            bcs = solver_session.settings.setup.boundary_conditions
            bcs.settings.physical_velocity_porous_formulation = True
            bcs.set_zone_type(zone_list=[self.INLET_NAME], new_type="velocity-inlet")
            inlet = bcs.velocity_inlet[self.INLET_NAME]
            inlet.momentum.velocity_magnitude                = self.setup.INLET_VELOCITY
            inlet.thermal.temperature.value                  = self.setup.INLET_TEMPERATURE
            inlet.species.species_mass_fraction["so2"].value = inlet_Y_SO2
            inlet.species.species_mass_fraction["o2"].value  = inlet_y_o2

            outlet = bcs.pressure_outlet[self.OUTLET_NAME]
            outlet.thermal.backflow_total_temperature.value = self.setup.OUTLET_BACKFLOW_TEMPERATURE

            rep_defs = ReportDefinitions(solver_session)

            for i in range(len(self.setup.RESPONSES)):

                rep_defs.surface.create(name=self.setup.RESPONSES[i])
                resp                  = rep_defs.surface[self.setup.RESPONSES[i]]
                resp.report_type      = self.setup.RESPONSES_SCOPE[i]
                resp.field            = self.setup.REPSONSES_FIELD[i]
                resp.surface_names    = [self.setup.RESPONSES_SURFACE[i]]

                solver_session.settings.solution.monitor.report_files.create(name=self.setup.RESPONSES[i])
                solver_session.settings.solution.monitor.report_files[self.setup.RESPONSES[i]] = {
                "file_name": str(pathlib.Path(self.cwd) / (self.setup.RESPONSES[i] + ".out")), "report_defs": [self.setup.RESPONSES[i]]}


            initialization = Initialization(solver_session)
            initialization.initialization_type = "hybrid"
            initialization.initialize()

            run_calc = RunCalculation(solver_session)
            run_calc.iter_count = self.setup.ITER_COUNT
            run_calc.calculate()
            solver_session.settings.file.write_case_data(file_name= case_path)
            print(f"  [3/3] Solver complete. Case saved to {case_path}")

            response_values = []
            for file in self.setup.RESPONSES:
                filepath = pathlib.Path(self.cwd) / (file + ".out")
                value = self.read_response(filepath= pathlib.Path(filepath))
                if value == False:
                    pass #Placeholder for convergence issue handling
                response_values.append(value)

            so3, so2 = response_values
            conversion = self._conversion(so3= so3, so2= so2)

            self.write_response([so2, conversion])
            return so2, conversion
        
        except Exception as e:
            print(f"  Solver failed: {type(e).__name__}: {e}")
            traceback.print_exc()
            raise
        finally:
            solver_session.exit()
        
        

    @staticmethod
    def _conversion(so3: float, so2: float) -> float:
        M_so3, M_so2 = 80.06, 64.07
        return (so3 / M_so3) / ((so3 / M_so3) + (so2 / M_so2))
    
    def read_response(self, filepath : pathlib.Path):
        if filepath.exists():
            lines = [l for l in filepath.read_text().splitlines() if l.strip()]
            converged = int(lines[-1].split()[0]) < self.setup.ITER_COUNT
            status    = "converged" if converged else "did NOT converge"
            print(f"Case in {filepath.parent} {status}!")
            if status == "converged":
                value = float(lines[-1].split()[1])
                return value
            else:
                return False
        else:ytrsaERTYUIO'
        
            print(f"Warning: Could not read {filepath}")
        return False
    def write_response(self, values: list):
        output_path = self.cwd / "output.csv"
        with open(output_path, 'w', newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(values)

class SimulationRun:
    def __init__(self,
                 geometry_setup    : GeometrySetup,
                 mesh_setup        : MeshSetup,
                 simulation_setup  : SimulationSetup,
                 cwd               : PathLike,
                 LengthCat1        : float,
                 LenghtCool        : float,
                 run_geometry      : bool = True,
                 run_meshing       : bool = True,
                 ) -> None:

        self.cwd              = pathlib.Path(cwd)
        self.geometry_setup   = geometry_setup
        self.mesh_setup       = mesh_setup
        self.simulation_setup = simulation_setup
        self.inlet_Y_SO2      = 0.07
        self.run_geometry     = run_geometry
        self.run_meshing      = run_meshing

        # Override catalyst lengths in geometry setup
        self.geometry_setup.lengths[1] = LengthCat1
        self.geometry_setup.lengths[2] = LenghtCool

        self.geometry_file = self.cwd / (geometry_setup.name + ".pmdb")
        self.mesh_file     = self.cwd / "geometry.msh.h5"

    def run(self) -> tuple:

        # --- Stage 1: Geometry ---
        if self.run_geometry:
            print("[1/3] Running geometry...")
            geo_script = GeometryScript(
                GeometrySetup = self.geometry_setup,
                cwd           = self.cwd,
            )
            self.geometry_file = geo_script._run_geometry()
            print(f"  Geometry saved to: {self.geometry_file}")
        else:
            print(f"[1/3] Skipping geometry. Using: {self.geometry_file}")
            if not self.geometry_file.exists():
                raise FileNotFoundError(f"Geometry file not found: {self.geometry_file}")

        # --- Stage 2: Meshing ---
        if self.run_meshing:
            print("[2/3] Running meshing...")
            mesh_script = MeshScript(
                GeometryScript = GeometryScript(self.geometry_setup, self.cwd),
                MeshSetup      = self.mesh_setup,
                geometry_file  = self.geometry_file,
            )
            self.mesh_file = mesh_script._run_meshing()
        else:
            print(f"[2/3] Skipping meshing. Using: {self.mesh_file}")
            if not self.mesh_file.exists():
                raise FileNotFoundError(f"Mesh file not found: {self.mesh_file}")

        # --- Stage 3: Simulation ---
        print("[3/3] Running simulation...")
        mesh_script_ref = MeshScript(
            GeometryScript = GeometryScript(self.geometry_setup, self.cwd),
            MeshSetup      = self.mesh_setup,
            geometry_file  = self.geometry_file,
        )
        sim_script = SimulationScript(
            mesh_file        = self.mesh_file,
            SimulationSetup  = self.simulation_setup,
            MeshScript       = mesh_script_ref,
            GeometrySetup    = self.geometry_setup
        )
        so3, so2 = sim_script._run_simulation()
        conversion = SimulationScript._conversion(so3=so3, so2=so2)

        print(f"\n  Results:")
        print(f"    SO3 outlet fraction : {so3:.6f}")
        print(f"    SO2 outlet fraction : {so2:.6f}")
        print(f"    SO2 conversion      : {conversion:.4%}")

        return so2, conversion

geometry_setup   = GeometrySetup()
mesh_setup       = MeshSetup()
simulation_setup = SimulationSetup()
base_cwd         = pathlib.Path.cwd()

def run_simulation(LengthCat1:  float,
                   LenghtCool:  float,
                   iteration:   int = 0,
                   particle:    int = 0) -> tuple:

    run_cwd = base_cwd / f"Iteration_{iteration:03d}" / f"Particle_{particle:03d}"
    run_cwd.mkdir(parents=True, exist_ok=True)

    run = SimulationRun(
        geometry_setup   = geometry_setup,
        mesh_setup       = mesh_setup,
        simulation_setup = simulation_setup,
        cwd              = run_cwd,
        LengthCat1       = LengthCat1,
        LenghtCool       = LenghtCool,
    )
    return run.run()

config    = PSOConfig(
    n_params           = 2,
    n_responses        = 2,
    x_lb               = [500, 500],
    x_ub               = [1000, 1000],
    simulation_factory = run_simulation,
)

optimizer = PSOOptimizer(config= config, initial_particles= np.random.uniform(size= (2,15)))
run = PSOOptimizer.run()