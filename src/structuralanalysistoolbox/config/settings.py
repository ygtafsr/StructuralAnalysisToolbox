
# src/structuralanalysistoolbox/config/settings.py
from importlib import resources
from pathlib import Path
import shutil

DB_FILENAME = "materialdatabase.json"

def ensure_material_database(target_dir: Path | None = None) -> Path:
    if target_dir is None:
        target_dir = Path.cwd()

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    target_path = target_dir / DB_FILENAME

    if not target_path.exists():
        with resources.files("structuralanalysistoolbox.materials").joinpath(DB_FILENAME).open("rb") as src:
            with open(target_path, "wb") as dst:
                shutil.copyfileobj(src, dst)

    return target_path


# Root of the installed package
PACKAGE_ROOT = Path(__file__).resolve().parents[1]

# Default material database path
DEFAULT_MAT_DATABASE_PATH = PACKAGE_ROOT / "materials" / "materialdatabase.json"

# Optional: helper to ensure the file exists
if not DEFAULT_MAT_DATABASE_PATH.exists():
    raise FileNotFoundError(f"Material database not found at: {DEFAULT_MAT_DATABASE_PATH}")
