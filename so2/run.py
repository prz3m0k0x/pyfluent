
from pathlib import Path
import subprocess, os, json

def extract_last_value(filepath):
    if os.path.exists(filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
                    
            valid_lines = [line for line in lines if line.strip()]
            if valid_lines:
                last_line = valid_lines[-1]  
                parts = last_line.split()
                return float(parts[1])
                        
            print(f"Warning: Could not read {filepath}")
            return 0.0
        
def conversion(so3, so2):
    M_so3 = 80.06
    M_so2 = 64.07
    return (so3 / M_so3) / ((so3 / M_so3) + (so2 / M_so2))


def write_response(filepath, outputs, names):
    response= {name : output for name, output in zip(names, outputs)}
    json_string = json.dumps(response)
    with open(filepath, "w") as f:
        json.dump(response, f)


def run(LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2, path, script_file, mixture_file, solver_file):
    py = r"C:\Users\praktyka\AppData\Local\Programs\Python\Python313\python.exe"
    script = r"C:\Users\praktyka\Desktop\Golasz\optislang\test2\test2.opr\solution.py"
    args = [str(LengthCat1), str(LengthCool), str(LengthCat2), str(inlet_Y_SO2), str(path), str(script_file), str(mixture_file), str(solver_file)]
    

    subprocess.run([py, script, *args], check=True, cwd= path)

    extracted_so3 = float(extract_last_value(os.path.join(path, "so3-outlet-fraction.out")))
    extracted_so2 = float(extract_last_value(os.path.join(path, "so2-outlet-fraction.out")))
    conv = conversion(so3= extracted_so3, so2= extracted_so2)
    responses = [extracted_so3, extracted_so2, conv]
    names = ['so3', 'so2', 'conversion']
    response= write_response(filepath=(os.path.join(path, "response.json")), outputs= responses, names= names)

    return response

path= Path(DESIGN_DIR)

script_file_name = "geometry_script_3d_python_v.scscript"
mixture_file_name = "reaction-mixture.scm"
solver_file_name = "solution.py"

opr_dir = None
for parent_level in [path, *path.parents]:
    # Look for any folder ending in .opr at this level
    opr_matches = list(parent_level.glob("*.opr"))
    
    if opr_matches and opr_matches[0].is_dir():
        opr_dir = opr_matches[0]
        break

if opr_dir:

    script_file = opr_dir / script_file_name
    mixture_file = opr_dir / mixture_file_name
    solver_file = opr_dir / solver_file_name
    
    print(f"Successfully found .opr dir: {opr_dir}")
else:
    print("Error: Could not locate the .opr directory in any parent folders.")

output = run(LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2, path, script_file, mixture_file, solver_file)