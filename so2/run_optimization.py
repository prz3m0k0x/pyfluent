import os
from ansys.optislang.core import Optislang, examples
from pathlib import Path

def print_node_info(node):
    name = node.get_name()
    type_ = node.get_type()
    status = node.get_status()
    print(name, type_, status)

def process_nodes(nodes):
    for node in nodes:
        print_node_info(node= node)
        if isinstance(node, System):
            process_nodes(node.get_nodes())

path = Path.cwd()
project_name = "so2-optimization.opf"
with Optislang(project_path= path / project_name) as osl:
    print(osl)
    example = examples.get_files("calculator_with_params")[1][0]
    print(example)
    project = osl.application.project

    root_system = project.root_system
    nodes = root_system.get_nodes()
    process_nodes(nodes)