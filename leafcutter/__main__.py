from pathlib import Path
from runpy import run_path

pkg_dir = Path(__file__).resolve().parent

def leafcutter_cluster():
    script_pth = pkg_dir / "clustering" / "leafcutter_cluster_regtools.py"
    run_path(str(script_pth), run_name="__main__")
    
def leafcutter_ds():
    script_pth = pkg_dir / "differential_splicing" / "leafcutter_ds.py"
    run_path(str(script_pth), run_name="__main__")