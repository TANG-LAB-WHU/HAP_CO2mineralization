"""Dependency checks must inspect artifacts, not just scheduler success."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('job_record',Path(__file__).resolve().parents[1]/'workflow/hpc/job_record.py')
jr=importlib.util.module_from_spec(spec);spec.loader.exec_module(jr)

class HPCProvenance(unittest.TestCase):
    def test_validated_dependency_rejects_changed_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);artifact=d/'data';artifact.write_bytes(b'explicit test fixture')
            manifest=d/'manifest.json';record={'status':'completed','validation':{'passed':True},'outputs':[{'path':'data','sha256':jr.sha(artifact)}]}
            manifest.write_text(json.dumps(record));jr.validated_upstream(manifest)
            artifact.write_bytes(b'changed fixture')
            with self.assertRaises(ValueError):jr.validated_upstream(manifest)

    def test_process_success_does_not_grant_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'manifest.json';p.write_text(json.dumps({'status':'completed','validation':{'passed':False}}))
            with self.assertRaises(ValueError):jr.validated_upstream(p)
            p.write_text(json.dumps({'status':'completed','validation':{'passed':True},'outputs':[]}))
            with self.assertRaises(ValueError):jr.validated_upstream(p)
