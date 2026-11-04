"""Failure-path tests: bad engine output must not become successful scientific evidence."""
import sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
from recovery_dynamics import thermo

HEADER='Step Time Temp PotEng KinEng TotEng Press\n'
ROWS='0 0 200 -10 1 -9 0\n10 0.0025 201 -10.005 1.005 -9 0\n'
END='Loop time of 0.1 on 1 procs for 10 steps with 2 atoms\nTotal wall time: 0:00:01\n'
class ThermoFailureTests(unittest.TestCase):
    def parse(self,text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.log';p.write_text(text,encoding='utf-8');return thermo(p)
    def test_complete_block_preserves_energy_and_time(self):
        x=self.parse(HEADER+ROWS+END);self.assertEqual(x.shape,(2,7));np.testing.assert_allclose(x[:,1],[0,.0025]);np.testing.assert_allclose(x[:,3]+x[:,4],x[:,5])
    def test_missing_program_termination_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS+END.split('Total wall time:')[0])
    def test_error_even_with_termination_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS+END+'ERROR: Lost atoms\n')
    def test_nonfinite_energy_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS.replace('-10.005','nan')+END)
    def test_malformed_row_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS.replace('1.005 -9 0','1.005 -9')+END)
    def test_duplicate_time_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS.replace('10 0.0025','10 0')+END)
    def test_missing_loop_end_rejected(self):
        with self.assertRaises(ValueError):self.parse(HEADER+ROWS+'Total wall time: 0:00:01\n')
if __name__=='__main__':unittest.main()
