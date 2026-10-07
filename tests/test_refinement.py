"""Contract tests use explicit toy fixtures, never fallback scientific datasets."""
import sys,tempfile,unittest
from pathlib import Path
import numpy as np
from ase import Atoms
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
import refine_pipeline as rp

class RefinementContracts(unittest.TestCase):
    def test_unsupported_sampling_plan_is_rejected(self):
        cfg={'sampling':{'replicas':8,'timestep_fs':.25,'steps_per_replica':16,'frame_stride':4,'discard_initial_steps':4}}
        rp.validate_sampling_profile(cfg)
        for key,value in [('replicas',9),('timestep_fs',.5),('discard_initial_steps',17)]:
            changed={'sampling':dict(cfg['sampling'],**{key:value})}
            with self.assertRaises(ValueError):rp.validate_sampling_profile(changed)

    def test_charge_gate_rejects_non_neutral_fragment(self):
        a=Atoms('OH',positions=[[0,0,0],[0,0,.96]],cell=[10,10,10],pbc=True)
        self.assertEqual(rp.chemistry(a)['formal_charge_e'],-1)
        self.assertFalse(rp.chemistry(a)['geometry_gate'])

    def test_minimum_contact_gate(self):
        a=Atoms('CO2',positions=[[0,0,0],[0,0,.2],[0,0,-1.2]],cell=[10,10,10],pbc=True)
        self.assertFalse(rp.chemistry(a)['geometry_gate'])

    def test_d3_is_sibling_of_functional(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);a=Atoms('CO2',positions=[[5,5,5],[5,5,6.2],[5,5,3.8]],cell=[10,10,10],pbc=True)
            rp.cp_input(a,d,'test',dispersion=True);stack=[];parents=[]
            for line in (d/'test.inp').read_text().splitlines():
                parts=line.strip().split()
                if not parts:continue
                if parts[0]=='&END':stack.pop()
                elif parts[0].startswith('&'):
                    if parts[0]=='&VDW_POTENTIAL':parents.append(stack[-1])
                    stack.append(parts[0])
            self.assertEqual(parents,['&XC']);self.assertEqual(stack,[])

    def test_actual_cp2k_force_format_and_unit_conversion(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'fixture.log'
            p.write_text('FORCES| Atomic forces [hartree/bohr]\nFORCES| Atom x y z |f|\nFORCES| 1 1.0E-2 -2.0E-2 3.0E-2 4.0E-2\n')
            f=rp.forces_from_log(p)
            np.testing.assert_allclose(f[0],np.array([[.01,-.02,.03]])*rp.HARTREE/rp.BOHR)

    def test_leakage_descriptor_ignores_rigid_translation(self):
        a=Atoms('CO2',positions=[[5,5,5],[5,5,6.2],[5,5,3.8]],cell=[15,15,15],pbc=True)
        b=a.copy();b.translate([1.7,-.3,.5])
        np.testing.assert_allclose(rp.pair_features(a),rp.pair_features(b),atol=1e-12)

if __name__=='__main__':unittest.main()
