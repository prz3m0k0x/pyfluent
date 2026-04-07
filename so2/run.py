from pathlib import Path
import subprocess
import os
import json
import sys


def extract_last_value(filepath):
    filepath = Path(filepath)
    if filepath.exists():
        with open(filepath, "r") as f:
            lines = [line for line in f.readlines() if line.strip()]
            if lines:
                parts = lines[-1].split()
                return float(parts[1])
    print(f"Warning: Could not read {filepath}")
    return 0.0


def conversion(so3, so2):
    M_so3 = 80.06
    M_so2 = 64.07
    return (so3 / M_so3) / ((so3 / M_so3) + (so2 / M_so2))


def write_response(filepath, outputs, names):
    response = {name: output for name, output in zip(names, outputs)}
    with open(filepath, "w") as f:
        json.dump(response, f, indent=2)
    return response


def find_opr_dir(start_path: Path):
    for parent_level in [start_path, *start_path.parents]:
        opr_matches = list(parent_level.glob("*.opr"))
        if opr_matches and opr_matches[0].is_dir():
            return opr_matches[0]
    return None


def run(LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2,
        path, script_file, mixture_file, solver_file):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    args = [
        str(LengthCat1),
        str(LengthCool),
        str(LengthCat2),
        str(inlet_Y_SO2),
        str(path),
        str(script_file),
        str(mixture_file),
        str(solver_file),
    ]

    result = subprocess.run(
        [sys.executable, str(solver_file), *args],
        check=False,
        capture_output=True,
        text=True,
        cwd=path,
    )

    if result.returncode != 0:
        print("=== solution.py failed ===")
        print(result.stdout)
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, result.args)

    extracted_so3 = extract_last_value(path / "so3-outlet-fraction.out")
    extracted_so2 = extract_last_value(path / "so2-outlet-fraction.out")
    conv = conversion(so3=extracted_so3, so2=extracted_so2)

    responses = [extracted_so3, extracted_so2, conv]
    names = ["so3", "so2", "conversion"]
    response = write_response(path / "response.json", responses, names)

    return response


if __name__ == "__main__":
    argv = sys.argv[1:]

    if len(argv) == 8:
        LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2 = map(float, argv[:4])
        path = Path(argv[4])
        script_file = Path(argv[5])
        mixture_file = Path(argv[6])
        solver_file = Path(argv[7])

    elif len(argv) == 4:
        LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2 = map(float, argv)

        design_dir = os.environ.get("DESIGN_DIR")
        if not design_dir:
            raise RuntimeError("DESIGN_DIR environment variable not set.")

        path = Path(design_dir)
        path.mkdir(parents=True, exist_ok=True)

        opr_dir = find_opr_dir(path)
        if opr_dir is None:
            raise RuntimeError(f"Could not locate the .opr directory from {path}")

        script_file = opr_dir / "geometry_script_3d_python_v.scscript"
        mixture_file = opr_dir / "reaction-mixture.scm"
        solver_file = opr_dir / "solution.py"

    else:
        raise RuntimeError(
            "Usage:\n"
            "python run.py L1 Lcool L2 Y_SO2 particle_dir script mixture solver\n"
            "or\n"
            "python run.py L1 Lcool L2 Y_SO2   with DESIGN_DIR set"
        )

    output = run(
        LengthCat1, LengthCool, LengthCat2, inlet_Y_SO2,
        path, script_file, mixture_file, solver_file
    )
    print(f"Run complete: {output}")