from __future__ import annotations

import runpy
import sys
from pathlib import Path

# Compatibility launcher for users running from `Examples/` (Windows/CMD, etc.)
REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "examples" / "run_example.py"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

runpy.run_path(str(SCRIPT), run_name="__main__")
