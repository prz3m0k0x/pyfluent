import ansys.fluent.core as pyfluent

from ansys.fluent.core import examples
from ansys.fluent.core import Precision, UIMode

import_file_name = examples.download_file("mixing_elbow.msh.h5", "pyfluent/mixing_elbow")

solver = pyfluent.launch_fluent(
    precision= Precision.DOUBLE,
    processor_count= 2,
    mode= "solver",
    ui_mode= UIMode.GUI,
    cleanup_on_exit= False
)

print(solver.get_fluent_version())