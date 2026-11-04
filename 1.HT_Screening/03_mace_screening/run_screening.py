from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"workflow"))
from local_compute import main
if __name__=="__main__":
    sys.argv.insert(1,"mace")
    main()
