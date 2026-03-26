from pathlib import Path
import subprocess, os

so3_output = 0.0
so2_output = 0.0

def extract_last_value(filepath):
    if os.path.exists(filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
                    
            valid_lines = [line for line in lines if line.strip()]
            if valid_lines:
                last_line = valid_lines[-1]
                       
                parts = last_line.split()
                if len(parts) >= 2:
                    return float(parts[1])
                        
            print(f"Warning: Could not read {filepath}")
            return 0.0

def run(LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2, path, script_file, mixture_file, solver_file):
    py = r"C:\Users\praktyka\AppData\Local\Programs\Python\Python313\python.exe"
    script = r"C:\Users\praktyka\Desktop\Golasz\PyFluent\pyfluent\so2\solution.py"
    args = [str(LengthCat1), str(LengthCool), str(LengthCat2), str(inlet_Y_SO2), str(path), str(script_file), str(mixture_file), str(solver_file)]
    

    subprocess.run([py, script, *args], check=True, cwd= path)

    extracted_so3 = extract_last_value(os.path.join(path, "so3-outlet-fraction.out"))
    extracted_so2 = extract_last_value(os.path.join(path, "so2-outlet-fraction.out"))

    return float(extracted_so3), float(extracted_so2)
 
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

raw_so3,raw_so2 = run(LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2, path, script_file, mixture_file, solver_file)

so3_output = float(raw_so3)
so2_output = float(raw_so2)