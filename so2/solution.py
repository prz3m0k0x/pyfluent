import ansys.geometry.core as geometry
import ansys.fluent.core as pyfluent
from ansys.fluent.core.solver import (
    Species,
    Viscous,
    Mesh,
    Energy,
    Initialization,
    RunCalculation,
    ReportDefinitions
)
import os, traceback

#OPTIMIZATION PARAMETERS
LengthCat1 = 800.0
LengthCool = 1500.0
LengthCat2 = 1600.0
inlet_Y_SO2 = 0.1

directory = os.path.dirname(__file__)

with os.scandir(directory) as d:
    for e in d:
        if e.is_file and e.name.endswith(".scscript"):
            script_path = e.path
            break

if script_path is None:
    print("Warning: No .scscript file found!")

script_args = {"LengthBed1"   : "500",   # mm
                "LengthCat1"  : str(LengthCat1),   # mm
                "LengthCool"  : str(LengthCool),   # mm
                "LengthCat2"  : str(LengthCat2),   # mm
                "LengthBed2"  : "2500",   # mm
                "Height"      : "2500",   # mm
                "dir"         : directory,
                "name"        : "geometry"
                }


def Geometry(file_path : str, script_args : dict):
    modeler = geometry.launch_modeler_with_spaceclaim()
    modeler.run_discovery_script_file(file_path= file_path, script_args= script_args)
    geometry_file = os.path.join(script_args.get("dir"), "geometry.pmdb")
    return modeler, geometry_file

def Mesher(
        filename: str, #.pmbd file of case
        directory: str, #path to the dir of case
        ):
    
    #MESH PARAMETERS
    GLOBAL_MIN: float = 60.
    GLOBAL_MAX: float = 240.
    LOCAL_SIZE: float = 120.


    filename_strip = filename[:-5]
    mesh_save_name = f"{filename_strip}.msh.h5"
    mesh_path = os.path.join(directory, mesh_save_name)


    meshing_session = pyfluent.launch_fluent(
        mode= "meshing",
        precision= "double",
        processor_count=8,
        ui_mode="no_gui",
        cwd= directory,
        cleanup_on_exit= True
    )

    try:
        workflow = meshing_session.workflow  
        workflow.InitializeWorkflow(WorkflowType = "Watertight Geometry")
        
        workflow.TaskObject["Import Geometry"].Arguments = dict(FileName = filename)
        workflow.TaskObject["Import Geometry"].Execute()

        local_sizing = meshing_session.workflow.TaskObject["Add Local Sizing"]
        local_sizing.Arguments.set_state(
            {
                "AddChild": "yes",
                "BOIControlName": "body_size1",
                "BOIFaceLabelList": ["fluid_catalyst2", "fluid_catalyst1"],
                "BOIExecution": "Body Size",
                "BOISize": LOCAL_SIZE,
            }
        )
        local_sizing.AddChildAndUpdate()

        surface_mesh_gen = meshing_session.workflow.TaskObject["Generate the Surface Mesh"]
        surface_mesh_gen.Arguments.set_state(
            {"CFDSurfaceMeshControls": {"MaxSize": GLOBAL_MAX, "MinSize": GLOBAL_MIN}}
        )

        surface_mesh_gen.Execute()        


        workflow.TaskObject["Describe Geometry"].Arguments.set_state({
            "SetupType": "The geometry consists of only fluid regions with no voids",
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
                "VolumeFill": "poly-hexcore",
                "VolumeMeshPreferences": {
                    "CheckSelfProximity": "yes",
                    "ShowVolumeMeshPreferences": True,
                },
            }
        )

        volume_mesh_gen.Execute()

        meshing_session.tui.file.write_mesh(mesh_save_name)

    except Exception as e:
        print(f"Meshing failed for {filename}")
        print("Type:", type(e).__name__)
        print("Repr:", repr(e))
        traceback.print_exc()
        return None

    finally:
        meshing_session.exit()
    
    return meshing_session, mesh_save_name

def Setup(
    filename: str,          #path to .msh.h5 file
    directory: str,     #working directory
    inlety_Y_SO2: float
):
    #CASE PARAMETERS

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

    VISCOUS_RESISTANCE  = 847936.476    #1/m^2
    INERTIAL_RESISTANCE = 283.104780 #1/m
    CATALYST_POROSITY = 0.6

    HEAT_EXCHANGER_VISCOUS_RESISTANCE = 1e4
    HEAT_EXCHANGER_INERTIAL_RESISTANCE = 100


    COOLING_HEAT_SINK   = -16000   #W/m^3

    # Inlet conditions
    INLET_VELOCITY = 0.5
    INLET_TEMPERATURE = 773.15 # K
    INLET_Y_SO2 = inlety_Y_SO2       #SO2 mass fraction at inlet
    INLET_Y_O2 = 0.21 * (1 - INLET_Y_SO2)   #O2 in remaining air

    # Solver
    ITER_COUNT          = 500

    solver_session = pyfluent.launch_fluent(
        mode=pyfluent.FluentMode.SOLVER,
        cwd=directory,
        dimension=3,
        precision=pyfluent.Precision.DOUBLE,
        processor_count=8,
        cleanup_on_exit=True,
        ui_mode="no_gui",
    )

    try:
        solver_session.settings.file.read_mesh(file_name=filename)

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

        general = solver_session.settings.setup.general
        general.solver.time.set_state("steady")
        general.solver.type.set_state("pressure-based")


        viscous = Viscous(solver_session)
        viscous.model = "k-epsilon"
        viscous.k_epsilon_model = "realizable"
        

        Energy(solver_session).enabled = True

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

        fluid_zone = solver_session.settings.setup.cell_zone_conditions.fluid
        for fluid in fluid_zone.get_object_names():
            fluid_zone[fluid].reaction.react = False


        # Porous zone setup for all catalytic porous zones
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
                'value': CATALYST_POROSITY}
            )
            fz.reaction.react = True

        # Volumetric heat sink in cooling zone
        fluid_zone[ZONE_COOLING].sources.set_state({"enable" : True})
        fluid_zone[ZONE_COOLING].sources.terms["energy"].resize(1)
        fluid_zone[ZONE_COOLING].sources.terms.set_state({'energy' : [{'option' : "value"}] })
        fluid_zone[ZONE_COOLING].sources.terms["energy"][0].set_state(COOLING_HEAT_SINK)
        fluid_zone[ZONE_COOLING].porous_zone.porous = True
        fluid_zone[ZONE_COOLING].porous_zone.viscous_resistance.set_state([
                    {"option": "value", "value": HEAT_EXCHANGER_VISCOUS_RESISTANCE},
                    {"option": "value", "value": HEAT_EXCHANGER_VISCOUS_RESISTANCE},
                    {"option": "value", "value": HEAT_EXCHANGER_VISCOUS_RESISTANCE}
        ])

        fluid_zone[ZONE_COOLING].porous_zone.inertial_resistance.set_state([
                    {"option": "value", "value": HEAT_EXCHANGER_INERTIAL_RESISTANCE},
                    {"option": "value", "value": HEAT_EXCHANGER_INERTIAL_RESISTANCE},
                    {"option": "value", "value": HEAT_EXCHANGER_INERTIAL_RESISTANCE}
        ])


        #Boundary conditions
        bcs = solver_session.settings.setup.boundary_conditions
        bcs.settings.physical_velocity_porous_formulation = True
        # Inlet
        bcs.set_zone_type(zone_list=[INLET_NAME], new_type="velocity-inlet")
        inlet = bcs.velocity_inlet[INLET_NAME]

        inlet.momentum.velocity_magnitude = INLET_VELOCITY
        inlet.thermal.temperature.value = INLET_TEMPERATURE
        inlet.species.species_mass_fraction["so2"].value = INLET_Y_SO2
        inlet.species.species_mass_fraction["o2"].value = INLET_Y_O2

        #outlet
        outlet = bcs.pressure_outlet[OUTLET_NAME]
        outlet.thermal.backflow_total_temperature.value = 700  # backflow temperature


        #Report definitions
        rep_defs = ReportDefinitions(solver_session)

        # SO3 mass flow rate at outlet
        rep_defs.surface.create(name="so3-outlet-fraction")
        so3_rep = rep_defs.surface["so3-outlet-fraction"]
        so3_rep.report_type = "surface-areaavg"
        so3_rep.field = 'so3'
        so3_rep.surface_names = [OUTLET_NAME]
        so3_rep.output_parameter = True

        # SO2
        rep_defs.surface.create(name="so2-outlet-fraction")
        so2_rep = rep_defs.surface["so2-outlet-fraction"]
        so2_rep.report_type = "surface-areaavg"
        so2_rep.field = "so2"
        so2_rep.surface_names = [OUTLET_NAME]

        #Initialization
        initialization = Initialization(solver_session)
        initialization.initialization_type = "hybrid"
        initialization.initialize()

        #Run
        run_calc = RunCalculation(solver_session)
        run_calc.iter_count = ITER_COUNT
        run_calc.calculate()

        so3_res = solver_session.settings.solution.report_definitions.compute(report_defs=["so3-outlet-fraction"])
        so2_res = solver_session.settings.solution.report_definitions.compute(report_defs=["so2-outlet-fraction"])

        extracted_so3 = so3_res[0]["so3-outlet-fraction"][0]
        extracted_so2 = so2_res[0]["so2-outlet-fraction"][0]

        # Save case and data
        case_path = os.path.join(directory, "so2_reactor.cas.h5")
        solver_session.settings.file.write_case_data(file_name=case_path)
        print(f"Setup complete. Case saved to {case_path}")

        return solver_session, extracted_so3, extracted_so2

    except Exception as e:
        print(f"Setup failed: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None

def Workflow(script_args, inlet_Y_SO2, directory, script_path):

    geometry = Geometry(file_path= script_path, script_args= script_args)[1]

    mesh = Mesher(filename= geometry, directory= directory)[1]

    solution = Setup(filename= "geometry.msh.h5", directory= directory, inlety_Y_SO2= inlet_Y_SO2)

    return solution[1], solution[2]

#Optimization responses

output = Workflow(script_args= script_args, inlet_Y_SO2= inlet_Y_SO2, directory= directory, script_path= script_path)
so3_mass_fraction_out = output[0]
so2_mass_fraction_out = output[1]

so3_conversion = so3_mass_fraction_out / inlet_Y_SO2

#Directory for geometry script

