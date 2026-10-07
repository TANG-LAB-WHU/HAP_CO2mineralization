import importlib.util, pathlib, unittest, tempfile, json, subprocess, sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

class LocalSafetyTests(unittest.TestCase):
    def test_no_small_dataset_leakage(self):
        split=load('split','3.ML_AIMD/Step3_mlip_deepmd/prepare_multisystem_input.py').split_systems_by_ratio
        for n in range(3):
            with self.assertRaises(ValueError): split([str(i) for i in range(n)],.8,.1,.1)
        for n in range(3,30):
            a,b,c=split([str(i) for i in range(n)],.8,.1,.1)
            self.assertEqual(len(set(a+b+c)),n)
            self.assertTrue(a and b and c)
            self.assertFalse(set(a)&set(b) or set(b)&set(c) or set(a)&set(c))
        with self.assertRaises(ValueError): split(['a','a','b'],.8,.1,.1)
        with self.assertRaises(ValueError): split(['a','b','c'],.8,-.1,.3)

    def test_missing_filter_inputs_fail(self):
        with tempfile.TemporaryDirectory() as d:
            r=subprocess.run([sys.executable,str(ROOT/'0.Generative_Discovery/post_filter.py'),'--results-dir',d],capture_output=True)
            self.assertNotEqual(r.returncode,0)

    def test_adsorption_reference_consistency(self):
        try: m=load('local_compute','workflow/local_compute.py')
        except ModuleNotFoundError: self.skipTest('numpy unavailable')
        refs={k:{'energy_eV':v,'model_sha256':'same'} for k,v in [('adsorbed',-14),('slab',-10),('co2',-3)]}
        self.assertEqual(m.adsorption_energy(refs),-1)
        refs['co2']['model_sha256']='different'
        with self.assertRaises(ValueError): m.adsorption_energy(refs)

    def test_cp2k_virial_sign_volume_units(self):
        import numpy as np
        m=load('cp2k_validate','workflow/cp2k_validate.py')
        stress=np.diag([1.,2.,3.])
        v=m.virial(stress,100.,'GPa')
        np.testing.assert_allclose(v,np.diag([.6241509074460763,1.2483018148921526,1.872452722338229]))
        np.testing.assert_allclose(v,m.virial(stress*10000,100,'bar'))
        # Tensile-positive ASE stress is negative CP2K stress: virial=-V*sigma_ASE.
        ase_stress=-stress*1e9/1.602176634e11
        np.testing.assert_allclose(v,-ase_stress*100)
        with self.assertRaises(ValueError):m.virial(stress,0,'bar')

    def test_incomplete_cp2k_rejected(self):
        m=load('cp2k_validate','workflow/cp2k_validate.py')
        with tempfile.TemporaryDirectory() as d:
            f=pathlib.Path(d)/'failed.log';f.write_text('SCF run NOT converged\nPROGRAM ENDED AT')
            with self.assertRaises(ValueError):m.parse(f)

    def test_legacy_cp2k_source_format(self):
        import ast,numpy as np
        path=ROOT/'3.ML_AIMD/Step2_dataset_synthesis/generate_deepmd_data.py'
        tree=ast.parse(path.read_text(encoding='utf-8'))
        func=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='parse_cp2k_virial_file')
        scope={'np':np};exec(compile(ast.Module(body=[func],type_ignores=[]),str(path),'exec'),scope)
        parse=scope['parse_cp2k_virial_file']
        with tempfile.TemporaryDirectory() as d:
            f=pathlib.Path(d)/'stress.dat'
            f.write_text('STRESS| Analytical stress tensor [bar]\nSTRESS| x y z\nSTRESS| x 1 0 0\nSTRESS| y 0 2 0\nSTRESS| z 0 0 3\nSTRESS| Trace 6\n')
            v=parse(f,1,np.array([np.eye(3)*10]))
            np.testing.assert_allclose(v.reshape(3,3),np.diag([1,2,3])*1000*1e5/1.602176634e11)
            with self.assertRaises(ValueError):parse(f,2,np.array([np.eye(3)*10]))
            f.write_text(f.read_text().replace('[bar]','[GPa]'))
            with self.assertRaises(ValueError):parse(f,1,np.array([np.eye(3)*10]))

if __name__=='__main__':unittest.main()
