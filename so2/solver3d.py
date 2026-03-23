import ansys.fluent.core as pyfluent
import os, traceback
from ansys.fluent.core.solver import (
    Species,
    Viscous,
    Mesh,
    Energy,
    Initialization,
    RunCalculation,
    ReportDefinitions,
    Monitor,
)

### CASE SETTINGS ###

INLET_NAME          = "inlet"
OUTLET_NAME         = "outlet"
WALLS_NAME          = "walls"
MIXTURE_SCM         = "reaction-mixture"
MIXTURE_MODEL       = "reaction-mixture"

ZONE_FLUID1         = "fluid_zone1"
ZONE_CAT1           = "fluid_catalyst1"
ZONE_COOLING        = "fluid_cooling"
ZONE_CAT2           = "fluid_catalyst2"
ZONE_FLUID2         = "fluid_zone2"

COOLING_ZONE        = [ZONE_COOLING]
CATALYST_ZONES      = [ZONE_CAT1, ZONE_CAT2]


VISCOUS_RESISTANCE  = 1e5    # 1/m^2
INERTIAL_RESISTANCE = 500   # 1/m


COOLING_HEAT_SINK   = -5e3   # W/m^3

# Inlet conditions
INLET_VELOCITY = 5
INLET_TEMPERATURE = 773.15     # K  (500°C)
INLET_Y_SO2 = 0.10       # SO2 mass fraction at inlet
INLET_Y_O2 = 0.21 * (1 - INLET_Y_SO2)   # O2 in remaining air



# Solver
ITER_COUNT          = 500


def Setup(
    mesh: str,          # path to .msh.h5 file
    directory: str,     # working directory
):
    solver_session = pyfluent.launch_fluent(
        mode=pyfluent.FluentMode.SOLVER,
        cwd=directory,
        dimension=3,
        precision=pyfluent.Precision.DOUBLE,
        processor_count=4,
        cleanup_on_exit=False,
        ui_mode="gui",
    )

    try:
        # ── Read mesh ─────────────────────────────────────────────────────────
        solver_session.settings.file.read_mesh(file_name=mesh)

        mesh_obj = Mesh(solver_session, new_instance_name="mesh")
        mesh_obj.surfaces_list = mesh_obj.surfaces_list.allowed_values()
        mesh_obj.display()

        # graphics = solver_session.settings.results.graphics
        # graphics.views.auto_scale()
        # if graphics.picture.use_window_resolution.is_active():
        #     graphics.picture.use_window_resolution = False
        # graphics.picture.x_resolution = 3840
        # graphics.picture.y_resolution = 2880
        # graphics.picture.save_picture(file_name="mesh.png")

        # ── General: steady, axisymmetric, pressure-based ─────────────────────
        general = solver_session.settings.setup.general
        general.solver.time.set_state("steady")
        general.solver.type.set_state("pressure-based")

        # ── Turbulence: laminar ───────────────────────────────────────────────
        viscous = Viscous(solver_session)
        viscous.model = "k-epsilon"
        viscous.k_epsilon_model = "realizable"
        

        # ── Energy: enabled (required for temperature-dependent reactions) ────
        Energy(solver_session).enabled = True

        # ── Species transport + reactions from .scm mixture definition ────────
        species = Species(solver_session)
        species.model.option = "species-transport"
        species.reactions.enable_volumetric_reactions = True

        scm_path = os.path.join(directory, f"{MIXTURE_SCM}.scm")
        try:
            solver_session.tui.define.materials.data_base.database_type(
                "user-defined",
                f"{scm_path}"
            )
            solver_session.tui.define.materials.copy("mixture", MIXTURE_MODEL)
            print(f"Loaded mixture from {scm_path}")
        except Exception as e:
            print(f"Failed to load mixture SCM: {e}")

        species.model.material = MIXTURE_MODEL

        # ── Cell zone conditions ──────────────────────────────────────────────
        fluid_zone = solver_session.settings.setup.cell_zone_conditions.fluid
        for fluid in fluid_zone.get_object_names():
            fluid_zone[fluid].reaction.react = False
            # if fluid in CATALYST_ZONES:
            #     continue
            # else:
            #     fluid_zone[fluid].reaction.react = False

        # Porous zone setup for all three porous zones
        for zone in CATALYST_ZONES:
            fz = fluid_zone[zone]
            fz.porous_zone.porous = True
            fz.porous_zone.viscous_resistance.set_state([
                    {"option": "value", "value": VISCOUS_RESISTANCE},
                    {"option": "value", "value": VISCOUS_RESISTANCE},
                    {"option": "value", "value": VISCOUS_RESISTANCE}
                ])
            fz.porous_zone.inertial_resistance.set_state([
                {"option": "value", "value": INERTIAL_RESISTANCE},
                {"option": "value", "value": INERTIAL_RESISTANCE},
                {"option": "value", "value": INERTIAL_RESISTANCE}
            ])
            # Porosity
            fz.porous_zone.porosity.set_state(
                {'option': 'value',
                'value': 0.4}
            )

        # Volumetric heat sink in cooling zone
        fluid_zone[ZONE_COOLING].sources.set_state({"enable" : True})
        fluid_zone[ZONE_COOLING].sources.terms["energy"].resize(1)
        fluid_zone[ZONE_COOLING].sources.terms.set_state({'energy' : [{'option' : "value"}] })
        fluid_zone[ZONE_COOLING].sources.terms["energy"][0].set_state(COOLING_HEAT_SINK)

        # ── Boundary conditions ───────────────────────────────────────────────
        bcs = solver_session.settings.setup.boundary_conditions

        # Inlet: pressure inlet
        bcs.set_zone_type(zone_list=[INLET_NAME], new_type="velocity-inlet")
        inlet = bcs.velocity_inlet[INLET_NAME]

        inlet.momentum.velocity_magnitude = INLET_VELOCITY
        inlet.thermal.temperature.value = INLET_TEMPERATURE
        inlet.species.species_mass_fraction["so2"].value = INLET_Y_SO2
        inlet.species.species_mass_fraction["o2"].value = INLET_Y_O2

        outlet = bcs.pressure_outlet[OUTLET_NAME]
        outlet.thermal.backflow_total_temperature.value = 700  # backflow temperature


        # ── Report definitions (optiSLang objectives) ─────────────────────────
        rep_defs = ReportDefinitions(solver_session)

        # SO3 mass flow rate at outlet
        rep_defs.surface.create(name="so3-outlet-massflow")
        so3_rep = rep_defs.surface["so3-outlet-massflow"]
        so3_rep.report_type = "surface-massflowrate"
        so3_rep.field = "so3"
        so3_rep.surface_names = [OUTLET_NAME]

        # Pressure drop: inlet total pressure - outlet total pressure
        rep_defs.surface.create(name="p-inlet-avg")
        p_in = rep_defs.surface["p-inlet-avg"]
        p_in.report_type = "surface-areaavg"
        p_in.field = "pressure"
        p_in.surface_names = [INLET_NAME]

        rep_defs.surface.create(name="p-outlet-avg")
        p_out = rep_defs.surface["p-outlet-avg"]
        p_out.report_type = "surface-areaavg"
        p_out.field = "pressure"
        p_out.surface_names = [OUTLET_NAME]

        # SO2 conversion monitor (bonus — useful for convergence check)
        rep_defs.surface.create(name="so2-outlet-massflow")
        so2_rep = rep_defs.surface["so2-outlet-massflow"]
        so2_rep.report_type = "surface-massflowrate"
        so2_rep.field = "so2"
        so2_rep.surface_names = [OUTLET_NAME]

        # ── Monitor report files ──────────────────────────────────────────────
        monitor = Monitor(solver_session)

        so3_monitor_path  = os.path.join(directory, "so3-outlet-massflow.out")
        dp_monitor_path   = os.path.join(directory, "pressure-monitors.out")

        monitor.report_files.create(name=so3_monitor_path)
        monitor.report_files[so3_monitor_path].report_defs = ["so3-outlet-massflow", "so2-outlet-massflow"]

        monitor.report_files.create(name=dp_monitor_path)
        monitor.report_files[dp_monitor_path].report_defs = ["p-inlet-avg", "p-outlet-avg"]

        # ── Initialization ────────────────────────────────────────────────────
        initialization = Initialization(solver_session)
        initialization.initialization_type = "hybrid"
        initialization.initialize()

        # ── Run ───────────────────────────────────────────────────────────────
        run_calc = RunCalculation(solver_session)
        run_calc.iter_count = ITER_COUNT
        run_calc.calculate()

        # ── Save case and data ────────────────────────────────────────────────
        case_path = os.path.join(directory, "so2_reactor.cas.h5")
        solver_session.settings.file.write_case_data(file_name=case_path)
        print(f"Setup complete. Case saved to {case_path}")

        return solver_session

    except Exception as e:
        print(f"Setup failed: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    directory = r"C:\Users\hp\Desktop\Me\Studia\CFD\so2"
    mesh_file = os.path.join(directory, "mixing_tank.msh.h5")
    Setup(mesh=mesh_file, directory=directory)
