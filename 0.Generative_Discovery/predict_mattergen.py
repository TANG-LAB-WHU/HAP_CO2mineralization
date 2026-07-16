#!/usr/bin/env python3  
"""  
Generative Discovery of Novel Apatite-Family Materials via MatterGen  
Explores arbitrary chemical systems (e.g., Ca-P-O-H-C, Ca-Sr-P-O-H) using  
diffusion-based conditional generation, with TRI2024-corrected stability evaluation.  
"""  
  
import os  
import sys  
from pathlib import Path  
import argparse  
import subprocess  
  
def main():  
    parser = argparse.ArgumentParser(  
        description='Predict Material Phase Space with MatterGen'  
    )  
    parser.add_argument(
        '--chemical-system',
        type=str,
        default='Ca-P-O-H-C',
        help='Chemical system to condition on (default: Ca-P-O-H-C)'
    )
    parser.add_argument(  
        '--output-path',  
        type=str,  
        default=None,  
        help='Output directory path (default: results/{chemical-system}/)'  
    )  
    parser.add_argument(
        '--structures-output-path',
        type=str,
        default=None,
        help='Path to save relaxed structures during evaluation (default: same as --output-path)'
    )
    parser.add_argument(  
        '--model-name',  
        type=str,  
        default='chemical_system',  
        choices=['chemical_system', 'chemical_system_energy_above_hull'],  
        help='Pre-trained model to use (default: chemical_system)'  
    )  
    parser.add_argument(  
        '--batch-size',  
        type=int,  
        default=64,  
        help='Number of structures to generate per batch (default: 64)'  
    )  
    parser.add_argument(  
        '--num-batches',  
        type=int,  
        default=50,  
        help='Number of batches to generate (default: 50)'  
    )  
    parser.add_argument(  
        '--guidance-factor',  
        type=float,  
        default=2.0,  
        help='Diffusion guidance factor (default: 2.0)'  
    )  
    parser.add_argument(  
        '--energy-above-hull',  
        type=float,  
        default=0.05,  
        help='Energy above hull threshold (only for chemical_system_energy_above_hull model, default: 0.05)'  
    )  
    parser.add_argument(  
        '--evaluate',  
        action='store_true',  
        help='Run evaluation after generation (default: False)'  
    )  
    parser.add_argument(  
        '--record-trajectories',  
        action='store_true',  
        default=True,  
        help='Record denoising trajectories (default: True)'  
    )  
    parser.add_argument(  
        '--repo-root',  
        type=str,  
        default='../mattergen',  
        help='Path to MatterGen repository root (default: ../mattergen)'  
    )  
    parser.add_argument(  
        '--mattersim-model',  
        type=str,  
        default='MatterSim-v1.0.0-5M.pth',  
        choices=['MatterSim-v1.0.0-1M.pth', 'MatterSim-v1.0.0-5M.pth'],  
        help='MatterSim model to use for relaxation (default: MatterSim-v1.0.0-5M.pth)'  
    )  
    parser.add_argument(  
        '--skip-relax',  
        action='store_true',  
        help='Skip structure relaxation during evaluation (requires pre-computed energies)'  
    )  
    parser.add_argument(  
        '--energies-path',  
        type=str,  
        default=None,  
        help='Path to pre-computed energies file (.npy format, only used with --skip-relax)'  
    )  
  
    args = parser.parse_args()  
  
    output_path = Path(args.output_path) if args.output_path else Path(f'results/{args.chemical_system}/')
    structures_output_path = (
        Path(args.structures_output_path)
        if args.structures_output_path
        else output_path
    )
    # Ensure we pass a file path (not a directory) to evaluator
    if structures_output_path.suffix.lower() == '.extxyz':
        structures_output_file = structures_output_path
    else:
        structures_output_file = structures_output_path / 'relaxed_structures.extxyz'

    # Calculate total structures
    total_structures = args.batch_size * args.num_batches

    # Print configuration
    print("=" * 60)
    print(f"{args.chemical_system} Material Phase Space Prediction")
    print("=" * 60)
    print(f"Model: {args.model_name}")
    print(f"Chemical System: {args.chemical_system}")
    print(f"Output path: {output_path}")
    print(f"Batch size: {args.batch_size}")
    print(f"Number of batches: {args.num_batches}")
    print(f"Total structures: {total_structures}")
    print(f"Diffusion steps: 1000 (fixed for D3PM models)")
    print(f"Guidance factor: {args.guidance_factor}")
    if args.evaluate:
        print(f"MatterSim model: {args.mattersim_model}")
        print(f"Skip relaxation: {args.skip_relax}")
        print(f"Structures output path: {structures_output_file}")
    print("=" * 60)
    print()

    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    # Ensure parent directory for the relaxed structures file exists
    structures_output_file.parent.mkdir(parents=True, exist_ok=True)
  
    # Prepare properties for conditioning  
    if args.model_name == 'chemical_system':  
        properties = f"{{'chemical_system': '{args.chemical_system}'}}"  
    elif args.model_name == 'chemical_system_energy_above_hull':  
        properties = f"{{'chemical_system': '{args.chemical_system}', 'energy_above_hull': {args.energy_above_hull}}}"  
  
    # Build generation command  
    generate_cmd = [  
        'mattergen-generate',  
        str(output_path),  
        f'--pretrained-name={args.model_name}',  
        f'--batch_size={args.batch_size}',  
        f'--num_batches={args.num_batches}',  
        f'--properties_to_condition_on={properties}',  
        f'--diffusion_guidance_factor={args.guidance_factor}',  
        f'--record_trajectories={str(args.record_trajectories)}'  
    ]  
  
    print("Executing generation command:")  
    print(' '.join(generate_cmd))  
    print()  
  
    # Execute generation  
    result = subprocess.run(generate_cmd, check=False)  
  
    if result.returncode != 0:  
        print(f"\nError: Generation failed with return code {result.returncode}")  
        sys.exit(1)  
  
    print("\n" + "=" * 60)  
    print("Generation complete!")  
    print("=" * 60)  
    print(f"Generated files located at: {output_path}")  
    print("  - generated_crystals_cif.zip: Crystal structures in CIF format")  
    print("  - generated_crystals.extxyz: Structure file in ExtXYZ format")  
    print("  - generated_trajectories.zip: Denoising trajectories")  

    # Quick statistics on generated structures
    structures_input = output_path / 'generated_crystals.extxyz'
    if structures_input.exists():
        try:
            from ase.io import read
            from collections import Counter
            frames = read(structures_input, index=':')
            print(f"  - Successfully loaded {len(frames)} generated structures")
            
            elements = Counter()
            for f in frames:
                elements.update(f.get_chemical_symbols())
            print(f"  - Element distribution: {dict(elements)}")
            if frames:
                sizes = [len(f) for f in frames]
                print(f"  - Atoms/structure: Min={min(sizes)}, Max={max(sizes)}, Avg={sum(sizes)/len(sizes):.1f}")
        except ImportError:
            print("  - (Install 'ase' to see detailed structural statistics)")
  
    # Run evaluation if requested  
    if args.evaluate:  
        print("\n" + "=" * 60)  
        print("Starting evaluation of generated structures...")  
        print("=" * 60)  
  
        # Use specified repository root  
        repo_root = Path(args.repo_root).resolve()  
  
        # Verify the path is valid  
        if not (repo_root / '.git').exists():  
            print(f"Error: {repo_root} is not a Git repository")  
            print("Please specify the correct path using --repo-root")  
            sys.exit(1)  
  
        print(f"Repository root: {repo_root}")  
  
        # Check if reference dataset exists  
        reference_dataset = repo_root / "data-release" / "alex-mp" / "reference_TRI2024correction.gz"  
  
        # Check if reference dataset exists and is a real file (not just an LFS pointer)  
        if not reference_dataset.exists() or reference_dataset.stat().st_size < 1024:  
            print("Checking reference dataset...")  
            print("Downloading reference dataset...")  
  
            # Change to repository root before running git lfs  
            original_dir = Path.cwd()  
            os.chdir(repo_root)  
  
            try:  
                subprocess.run([  
                    'git', 'lfs', 'pull', '-I',  
                    'data-release/alex-mp/reference_TRI2024correction.gz',  
                    '--exclude='  
                ], check=True)  
                
                # Verify if the pull actually succeeded
                if not reference_dataset.exists() or reference_dataset.stat().st_size < 1024:
                    print("Error: Dataset is missing or still looks like an LFS pointer after pull.")
                    print("Check your Git LFS installation and repository access.")
                    sys.exit(1)
                    
                print("Reference dataset downloaded successfully")  
            except subprocess.CalledProcessError as e:  
                print(f"Error downloading reference dataset: {e}")  
                print("Evaluation cannot proceed without reference dataset")  
                sys.exit(1)  
            finally:  
                # Change back to original directory  
                os.chdir(original_dir)  
  
        # Build evaluation command based on skip-relax option  
        evaluate_cmd = [
            'mattergen-evaluate',  
            f'--structures_path={structures_input}',
            f'--relax={not args.skip_relax}',  
            '--structure_matcher=disordered',  
            f'--save_as={output_path}/metrics.json',
            f'--save_detailed_as={output_path}/detailed_metrics.json',
            f'--structures_output_path={structures_output_file}',
            f'--reference_dataset_path={reference_dataset}',
            '--energy_correction_scheme=TRI2024'
        ]  
  
        # Add energies path if skip-relax is enabled and energies are provided  
        if args.skip_relax and args.energies_path:  
            evaluate_cmd.append(f'--energies_path={args.energies_path}')  
  
        # Add MatterSim model selection if not skipping relaxation  
        if not args.skip_relax:  
            evaluate_cmd.append(f'--potential_load_path={args.mattersim_model}')  
  
        print("\nExecuting evaluation command:")  
        print(' '.join(evaluate_cmd))  
        print()  
  
        result = subprocess.run(evaluate_cmd, check=False)  
  
        if result.returncode == 0:  
            print("\nEvaluation complete!")  
            print(f"Evaluation results saved at: {output_path}/metrics.json")  
        else:  
            print(f"\nWarning: Evaluation failed with return code {result.returncode}")  
  
    print("\n" + "=" * 60)  
    print("Script execution complete!")  
    print("=" * 60)  
  
if __name__ == '__main__':  
    main()