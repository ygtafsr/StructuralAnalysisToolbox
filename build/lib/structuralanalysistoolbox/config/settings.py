
from pathlib import Path
import os

# Root of the installed package
PACKAGE_ROOT = Path(__file__).resolve().parents[1]

# Default material database path
DEFAULT_MAT_DATABASE_PATH = PACKAGE_ROOT / "materials" / "materialdatabase.json"

# Optional: helper to ensure the file exists
if not DEFAULT_MAT_DATABASE_PATH.exists():
    raise FileNotFoundError(f"Material database not found at: {DEFAULT_MAT_DATABASE_PATH}")


