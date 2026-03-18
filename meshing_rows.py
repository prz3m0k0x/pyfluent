import os
import traceback
import ansys.fluent.core as pyfluent


def Mesher(
    filename: str,
    directory: str,
    global_min: float = 0.3,
    global_max: float = 1.2,
    local_size: float = 0.15,
):
    geometry_file = filename
    save_dir = rf"C:\Users\hp\Desktop\Me\Studia\CFD\adsorpcja2026\{directory}"
    mesh_save_name = "mesh2d.msh.h5"
    mesh_path = os.path.join(save_dir, mesh_save_name)

    os.makedirs(save_dir, exist_ok=True)

    meshing_session = pyfluent.launch_fluent(
        mode=pyfluent.FluentMode.MESHING,
        precision=pyfluent.Precision.DOUBLE,
        processor_count=4,
        ui_mode="gui",
    )

    try:  
        xy_meshing = meshing_session.two_dimensional_meshing()

        load_cad = xy_meshing.load_cad_geometry_2d
        load_cad.file_name = geometry_file
        load_cad.length_unit = "mm"
        load_cad()

        update_boundaries = xy_meshing.update_boundaries_2d
        update_boundaries.selection_type = "zone"
        update_boundaries()

        global_sizing = xy_meshing.define_global_sizing_2d
        global_sizing.max_size = global_max
        global_sizing.min_size = global_min
        global_sizing()

        add_local_sizing = xy_meshing.add_local_sizing_2d
        add_local_sizing.add_child = "No"
        add_local_sizing.boi_control_name = "edgesize_1"
        add_local_sizing.boi_execution = "Edge Size"
        add_local_sizing.boi_size = local_size
        add_local_sizing.boi_zoneor_label = "label"
        add_local_sizing.draw_size_control = True
        add_local_sizing.edge_label_list = ["adsorbent_walls"]
        add_local_sizing.add_child_and_update(defer_update=False)

        generate_surface_mesh = xy_meshing.generate_initial_surface_mesh
        mesh_preferences = generate_surface_mesh.surface_2d_preferences
        mesh_preferences.show_advanced_options = True
        mesh_preferences.merge_edge_zones_based_on_labels = "yes"
        mesh_preferences.merge_face_zones_based_on_labels = "yes"
        generate_surface_mesh()

        tasks = meshing_session.workflow.TaskObject
        export_mesh = tasks["Export Fluent 2D Mesh"]
        export_mesh.Arguments.set_state(
            {
                "FileName": mesh_path,
            }
        )
        export_mesh.Execute()

        print(f"Meshing {geometry_file} successful -> {mesh_path}")
        return mesh_path

    except Exception as e:
        print(f"Meshing failed for {geometry_file}")
        print("Type:", type(e).__name__)
        print("Repr:", repr(e))
        traceback.print_exc()
        return None

    finally:
        meshing_session.exit(wait=True)


n = 6
geometry_list = [f'geometria-{2**(i+1)}rows.pmdb' for i in range(n)]
directory = [f'geometria-{2**(i+1)}' for i in range(n)]
print(geometry_list)
print(directory)
for i,j in zip(geometry_list, directory):
    try:
        print(f'Processing file: {i}')
        Mesher(filename= i, directory= j)
        print(f"File {i} processed!")
    except:
        print(f"File {i} failed to process!")

# Mesher(filename= geometry_list[5], directory= directory[5])
print("Mehsing finnished")
