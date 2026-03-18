import ansys.fluent.core as pyfluent
import os
import re
import sys
from io import StringIO

UDF_IGNORE = ["adsorpcja_byname.c"]
WALL_ZONE = "adsorbent_walls"
FLUID_ZONE = "geometria-powierzchnia"
INLET_NAME = "inlet"
OUTLET_NAME = "outlet"
MIXTURE_MODEL = "water_rhodamine.scm"
ADSORBED_SPECIE_FORMULA = 'rh-b'
UDF_FOLDER = "UDF"
UDF_LIBRARY_NAME = "libudf"
EXECUTE_ON_DEMAND = ['clear_memo', 'facethread_memo', 'init_adsorption_cache']
EXECUTE_AT_THE_END = ['update_adsorption']
SOURCES = ["mass_source"]
TIME_STEP_COUNT = 1
ITER_COUNT = 5
TIME_STEP_SIZE = 0.05

from ansys.fluent.core.solver import (
    Species,
    Viscous,
    Mesh,
    Energy,
    VelocityInlet,
    Initialization,
    RunCalculation,
    ReportDefinitions,
    Monitor
)


def update_udf_zone_ids(udf_file_path, wall_zone_id, fluid_zone_id):

    
    # Read the entire file
    with open(udf_file_path, 'r') as f:
        content = f.read()
    
    # Replace the #define lines
    try:
        content = re.sub(
            r'#define\s+WALL_ZONE_ID\s+\d+',
            f'#define WALL_ZONE_ID {wall_zone_id}',
            content
        )
        content = re.sub(
            r'#define\s+FLUID_ZONE_ID\s+\d+',
            f'#define FLUID_ZONE_ID {fluid_zone_id}',
            content
        )
    
    
        # Write back to the same file
        with open(udf_file_path, 'w') as f:
            f.write(content)
    
        print(f"Updated {udf_file_path}: WALL_ZONE_ID={wall_zone_id}, FLUID_ZONE_ID={fluid_zone_id}")

    except:
        print(f"No fluid and wallzone dependencies in udf {udf_file_path}")



def ConcentrationToMassFrac(concentration, density, molar_mass):
    return concentration * molar_mass / density

def SolverSettings(mesh_file: str):
    #Initializing case
    file_path = os.path.dirname(os.path.realpath(__file__))
    print(file_path)

    solver_session = pyfluent.launch_fluent(

        mode=pyfluent.FluentMode.SOLVER,
        cwd= file_path,
        dimension=2,
        precision=pyfluent.Precision.DOUBLE,
        processor_count=4,
        ui_mode="gui",
        cleanup_on_exit=False,
    )
    #Reading mesh
    solver_session.settings.file.read_mesh(file_name = mesh_file)
    mesh = Mesh(solver_session, new_instance_name="mesh")
    mesh.surfaces_list = mesh.surfaces_list.allowed_values()
    mesh.display()
    #saving mesh picture
    graphics = solver_session.settings.results.graphics
    graphics.views.auto_scale()
    if graphics.picture.use_window_resolution.is_active():

        graphics.picture.use_window_resolution = False
 
    graphics.picture.x_resolution = 3840
    graphics.picture.y_resolution = 2880
    graphics.picture.save_picture(file_name="mesh.png")

    #general settings
    time = solver_session.settings.setup.general.solver.time
    time.set_state("unsteady-2nd-order")

    #turbulence model
    viscous = Viscous(solver_session)
    viscous.model = "laminar"

    #importing materials
    species = Species(solver_session)
    species.model.option = "species-transport"

    database = MIXTURE_MODEL
    scm_file = str(os.path.abspath(database))
    print("SCM path:", scm_file)
    print("SCM exists:", os.path.exists(scm_file))
 
    try:

        solver_session.tui.define.materials.data_base.database_type(
            "user-defined",
            f'{database}'
        )

        print("User-defined database loaded.")
 
        solver_session.tui.define.materials.copy("mixture", "my-mixture")
    except Exception as e:
        print(f"Failure! Material database not loaded: {e}")


    #Setting new mixture in cell zone
    
    fluid_zone = solver_session.settings.setup.cell_zone_conditions.fluid
    fluid_zone_list = fluid_zone.get_object_names()
    fluid_zone_name = fluid_zone_list[0]
    fluid_zone.rename(FLUID_ZONE, f'{fluid_zone_name}')

    fluid_zone_name = FLUID_ZONE

    try:
        species.model.material = "my-mixture"
    except:
        print("Couldn't set mixture")

    #turning off energy as it is not neccesary
    Energy(solver_session).enabled = False

    #setting inlet conditions
    solver_session.settings.setup.boundary_conditions.set_zone_type(zone_list=[INLET_NAME], new_type= "velocity-inlet")
    inlet = VelocityInlet(solver_session, name="inlet")

    inlet.momentum.velocity_magnitude = 0.01
    inlet.species.species_mass_fraction[ADSORBED_SPECIE_FORMULA] = ConcentrationToMassFrac(0.005, 998.2, 470.0)

    wall_id = solver_session.scheme_eval.scheme_eval(f'(thread-name->id "{WALL_ZONE}")')
    fluid_id = solver_session.scheme_eval.scheme_eval(f'(thread-name->id "{FLUID_ZONE}")')

    user_defined = solver_session.settings.setup.user_defined
    user_defined.memory.memory_locations.set_state(5)

    # directory of this Python file
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # UDF folder inside it
    udf_dir = os.path.join(base_dir, UDF_FOLDER)
    
    c_files = []
    for root, dirs, files in os.walk(udf_dir):
        for name in files:
            if name.endswith(".c"):
                if name in UDF_IGNORE:
                    continue
                c_files.append(os.path.join(root, name))
    for file in c_files:
        path = os.path.join(udf_dir, file)
        update_udf_zone_ids(path, wall_zone_id= wall_id, fluid_zone_id= fluid_id)
    
    print(c_files)

    library_path = os.path.join(file_path, UDF_LIBRARY_NAME)


    user_defined.compiled_udf(library_name= library_path, source_files= c_files, use_built_in_compiler= True)
    user_defined.load(udf_library_name= library_path)

    user_defined.function_hooks.execute_at_end(lib_name= f"{EXECUTE_AT_THE_END[0]}::{UDF_LIBRARY_NAME}")

    fluid_zone[fluid_zone_name].sources.enable = True
    fluid_zone[fluid_zone_name].sources.terms["mass"].resize(1)
    fluid_zone[fluid_zone_name].sources.terms["species-0"].resize(1)

    fluid_zone[fluid_zone_name].sources.terms["mass"][0].option.set_state('udf')
    fluid_zone[fluid_zone_name].sources.terms['mass'][0].udf.set_state(f"mass_source::{UDF_LIBRARY_NAME}")

    fluid_zone[fluid_zone_name].sources.terms["species-0"][0].option.set_state('udf')
    fluid_zone[fluid_zone_name].sources.terms["species-0"][0].udf.set_state(f"mass_source::{UDF_LIBRARY_NAME}")

    initialization = Initialization(solver_session)

    initialization.initialization_type = "hybrid"
    initialization.initialize()

    

    for on_demand in EXECUTE_ON_DEMAND:
        user_defined.execute_on_demand(lib_name= f"{on_demand}::{UDF_LIBRARY_NAME}")

        rep_defs = ReportDefinitions(solver_session)
    
    rep_defs.surface.create(name= "specie0-mass-frac")
    surf_rep = rep_defs.surface["specie0-mass-frac"]
    surf_rep.report_type = "surface-areaavg"
    surf_rep.field = f"{ADSORBED_SPECIE_FORMULA}"
    surf_rep.surface_names = ["outlet"]

 
    rep_defs.volume.create(name= "adsorption-in-bed")
    volume_rep = rep_defs.volume["adsorption-in-bed"]
    volume_rep.report_type = "volume-massintegral"
    volume_rep.field = "udm-4"
    volume_rep.cell_zones = [f"{fluid_zone_name}"] 

    monitor = Monitor(solver_session)
    species_path = os.path.join(file_path, "specie0-mass-frac")
    adsorption_path = os.path.join(file_path, "adsorption-in-bed")

    monitor.report_files.create(name= species_path)
    monitor.report_files[species_path].report_defs = ["specie0-mass-frac"]

    monitor.report_files.create(name= adsorption_path)
    monitor.report_files[adsorption_path].report_defs = ["adsorption-in-bed"]

    transient_controls = solver_session.settings.solution.run_calculation.transient_controls
    
    transient_controls.time_step_count = TIME_STEP_COUNT
    transient_controls.time_step_size = TIME_STEP_SIZE
    transient_controls.max_iter_per_time_step = ITER_COUNT
    
    run_calc = RunCalculation(solver_session)
    run_calc.calculate()


    # Capture the console output
    old_stdout = sys.stdout
    sys.stdout = captured = StringIO()

    solver_session.settings.parallel.timer.usage()  # This prints your stats

    sys.stdout = old_stdout  # Restore console
    timer_stats = captured.getvalue()  # Get the captured text

    print("Captured:", timer_stats[:100] + "...")  # Verify it worked

    with open(os.path.join(file_path, "final_performance.txt"), "w") as f:
        f.write(timer_stats)

    print("Saved to final_performance.txt")

    return solver_session
 
 
if __name__ == "__main__":

    solver = SolverSettings(mesh_file="mesh2d.msh.h5")

