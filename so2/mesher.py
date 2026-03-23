import ansys.fluent.core as pyfluent
import os, traceback

def Mesher(
    filename: str, #.pmbd file of case
    directory: str, #path to the dir of case
    global_min: float = 1,
    global_max: float = 10,
    # local_size: float = 1,
    ):
    filename_strip = filename[:-5]
    mesh_save_name = f"{filename_strip}.msh.h5"
    mesh_path = os.path.join(directory, mesh_save_name)


    meshing_session = pyfluent.launch_fluent(
        mode= "meshing",
        precision= "double",
        processor_count=4,
        ui_mode="gui",
        cwd= directory,
        cleanup_on_exit= False
    )

    try:
        workflow = meshing_session.workflow  
        workflow.InitializeWorkflow(WorkflowType = "Watertight Geometry")
        
        workflow.TaskObject["Import Geometry"].Arguments = dict(FileName = filename)
        workflow.TaskObject["Import Geometry"].Execute()

        workflow.TaskObject["Describe Geometry"].Arguments.set_state({
            "SetupType": "The geometry consists of only fluid regions with no voids",
            "WallToInternal": "Yes",
        })
        workflow.TaskObject["Describe Geometry"].Execute()


        workflow.TaskObject["Update Boundaries"].Execute()

        workflow.TaskObject["Update Regions"].Execute()

        workflow.TaskObject["Add Boundary Layers"].InsertCompoundChildTask()
        workflow.TaskObject["smooth-transition_1"].Execute()

        workflow.TaskObject["Generate the Volume Mesh"].Execute()

        meshing_session.tui.file.write_mesh("mixing_tank.msh.h5")


    #     load_cad = xy_meshing.load_cad_geometry_2d
    #     load_cad.file_name = filename
    #     load_cad.refacet = "no"
    #     load_cad.length_unit = "mm"
    #     load_cad()

    #     update_boundaries = xy_meshing.update_boundaries_2d
    #     update_boundaries.selection_type = "label"
    #     update_boundaries()

    #     global_sizing = xy_meshing.define_global_sizing_2d
    #     global_sizing.max_size = global_max
    #     global_sizing.min_size = global_min
    #     global_sizing()

    #     # add_local_sizing = xy_meshing.add_local_sizing_2d
    #     # add_local_sizing.add_child = "yes"
    #     # add_local_sizing.boi_control_name = "facesize_1"
    #     # add_local_sizing.boi_execution = "Face Size"
    #     # add_local_sizing.boi_zoneor_label = "label"
    #     # add_local_sizing.boi_face_label_list = ["fluid_catalyst1", "fluid_catalyst2"]
    #     # add_local_sizing.boi_size = local_size
    #     # add_local_sizing.draw_size_control = True
    #     # add_local_sizing.add_child_and_update(defer_update=False)

    #     boundary_layers = xy_meshing.add_2d_boundary_layers
    #     boundary_layers.add_child = "yes"
    #     boundary_layers.control_name = "walls-inflation"
    #     boundary_layers.number_of_layers = 5
    #     boundary_layers.offset_method_type = "uniform"
    #     boundary_layers.add_child_and_update(defer_update=False)
    #     boundary_layers.execute()

    #     generate_surface_mesh = xy_meshing.generate_initial_surface_mesh

    #     mesh_preferences = generate_surface_mesh.surface_2d_preferences
    #     generate_surface_mesh.generate_quads = True
    #     mesh_preferences.show_advanced_options = True
    #     mesh_preferences.merge_edge_zones_based_on_labels = "yes"
    #     mesh_preferences.merge_face_zones_based_on_labels = "yes"
    #     generate_surface_mesh()

    #     tasks = meshing_session.workflow.TaskObject
    #     export_mesh = tasks["Export Fluent 2D Mesh"]
    #     export_mesh.Arguments.set_state(
    #         {
    #             "FileName": mesh_path,
    #         }
    #     )
    #     export_mesh.Execute()

    #     print(f"Meshing {filename} successful -> {mesh_path}")
    #     return mesh_path

    except Exception as e:
        print(f"Meshing failed for {filename}")
        print("Type:", type(e).__name__)
        print("Repr:", repr(e))
        traceback.print_exc()
        return None

    finally:
        meshing_session.exit()

Mesher("geometria.pmdb", r"C:\Users\hp\Desktop\Me\Studia\CFD\so2")