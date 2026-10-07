from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from screening_tools import main
if __name__ == "__main__": main('energy')
