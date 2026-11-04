"""One-sample official MatterGen smoke test with explicit resource/provenance limits."""
import argparse,hashlib,json,os,shutil,time
from pathlib import Path
import torch
p=argparse.ArgumentParser();p.add_argument('stage',choices=['generate','evaluate']);p.add_argument('--runs',type=Path,default=Path('.research-runs'));a=p.parse_args()
root=a.runs.resolve(); out=root/'mattergen-smoke';out.mkdir(exist_ok=True)
manifest_path=out/f'{a.stage}.manifest.json'
if manifest_path.exists():
    previous=json.loads(manifest_path.read_text())
    if previous.get('status')=='completed':raise FileExistsError('Preserving completed stage manifest')
    shutil.copyfile(manifest_path,out/f'{a.stage}-previous-{time.time_ns()}.manifest.json')
torch.set_num_threads(4)
if not torch.cuda.is_available():raise RuntimeError('CUDA required')
torch.cuda.set_per_process_memory_fraction(.35)
start=time.time()
record={'stage':a.stage,'status':'running','seed':42,'batch_size':1,'num_batches':1,'chemical_system':'Ca-P-O-H-C','manuscript_eligible':False}
record['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
manifest_path.write_text(json.dumps(record,indent=2))
def check(path,expected):
    actual=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
    if actual!=expected:raise ValueError('Input checksum mismatch: '+str(path))
    return actual
try:
    if a.stage=='generate':
        if (out/'generated_crystals.extxyz').exists():raise FileExistsError('Preserving generated sample')
        record['checkpoint_sha256']=check(root/'mattergen/checkpoints/chemical_system/checkpoints/last.ckpt','4ad21d977c2b776a31c92498c4d1c88d3e1e286deeb941e997af8e082310b80e')
        from mattergen.scripts.generate import main
        main(str(out),model_path=str(root/'mattergen/checkpoints/chemical_system'),batch_size=1,num_batches=1,properties_to_condition_on={'chemical_system':'Ca-P-O-H-C'},diffusion_guidance_factor=2.0,record_trajectories=False,seed=42)
    else:
        # The pinned reference encodes the newer pymatgen module path for the same
        # ComputedEntry classes. MatterGen's monty pin resolves pymatgen 2024, where
        # those classes live in entries.computed_entries. Decode through that alias.
        import sys
        import pymatgen.entries.computed_entries as computed_entries
        sys.modules.setdefault('pymatgen.core.entries',computed_entries)
        from mattergen.evaluation.reference import reference_dataset_serializer as serializer
        # Upstream decompress reads the entire uncompressed DB into RAM. Stream identical bytes
        # to disk instead, avoiding a multi-GB transient allocation. No scientific data changed.
        import gzip
        def streaming_decompress(path,directory):
            dest=Path(directory)/Path(path).name[:-3]
            with gzip.open(path,'rb') as src, dest.open('wb') as dst:shutil.copyfileobj(src,dst,16*1024*1024)
            return dest
        serializer.gzip_decompress=streaming_decompress
        # LMDB 3 rejects reopening the same environment in one process. Upstream's
        # index builder does so for each metadata key. Read identical keys through
        # its already-open read-only transaction instead; preserve all entries.
        def index_using_open_environment(self,path):
            with self.env.begin() as txn:
                return {system:{formula:serializer.lmdb_get(txn,f'{system}.{formula}.length')
                    for formula in serializer.lmdb_get(txn,f'{system}.reduced_formulas')}
                    for system in serializer.lmdb_get(txn,'chemical_systems')}
        serializer.LMDBBackedReferenceDatasetImpl._build_num_entries_by_chemsys_reduced_formulas=index_using_open_environment
        record['compatibility_fixes']=['streaming lossless gzip decompression','reuse read-only LMDB environment for metadata index','pymatgen.core.entries serialization alias to installed computed_entries']
        reference=root/'mattergen/data-release/alex-mp/reference_TRI2024correction.gz'
        record['reference_sha256']=check(reference,'ab1b4426901ff7b0f8af6ee9c3cc0dd2fbe4b65c758ef2a95683071492d1c51d')
        potential=root/'models/mattersim-v1.0.0-1M.pth'
        record['mattersim_sha256']=check(potential,'28b0b0b0f13efefee06b47ea4c9105a26bd3e2c8396da193430da96b3b49a8be')
        from mattergen.scripts.evaluate import main
        main(str(out/'generated_crystals.extxyz'),relax=True,save_as=str(out/'metrics.json'),save_detailed_as=str(out/'detailed_metrics.json'),potential_load_path=str(potential),reference_dataset_path=str(reference),device='cuda',structures_output_path=str(out/'relaxed_structures.extxyz'),energy_correction_scheme='TRI2024')
    record['status']='completed'
except Exception as e:
    record['status']='failed';record['error']=repr(e);raise
finally:
    record['elapsed_seconds']=time.time()-start
    record['peak_cuda_allocated_bytes']=torch.cuda.max_memory_allocated()
    (out/f'{a.stage}.manifest.json').write_text(json.dumps(record,indent=2))
