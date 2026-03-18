import ansys.fluent.core as pyfluent

import os
 
from ansys.fluent.core.solver import (
    Species,
    Viscous,
    Mesh,
    Energy,
    VelocityInlet,
    FluidCellZones
)

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

    database = "water_rhodamine.scm"
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
    fluid_zone.rename('geometria-powierzchnia', f'{fluid_zone_name}')

    fluid_zone_name = fluid_zone_name

    try:
        species.model.material = "my-mixture"
    except:
        print("Couldn't set mixture")

    #turning off energy as it is not neccesary
    Energy(solver_session).enabled = False

    #setting inlet conditions
    solver_session.settings.setup.boundary_conditions.set_zone_type(zone_list=["inlet"], new_type= "velocity-inlet")
    inlet = VelocityInlet(solver_session, name="inlet")

    inlet.momentum.velocity_magnitude = 0.01
    inlet.species.species_mass_fraction["rh-b"] = ConcentrationToMassFrac(0.005, 998.2, 470.0)

    user_defined = solver_session.settings.setup.user_defined
    user_defined.memory.memory_locations.set_state(5)
    UDF_names = ["test_UDF.c"]

    user_defined.compiled_udf(library_name= "libudf", source_files= UDF_names, use_built_in_compiler= True)


    return solver_session
 
 
if __name__ == "__main__":

    solver = SolverSettings(mesh_file="mesh2d.msh.h5")
 