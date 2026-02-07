#!/usr/bin/env python3
"""
Dynamic configuration generator for HAP workflow
Generates configuration files based on HKL parameter from folder name
Maintains English-only text as per user preference
"""

import os
import sys
import shutil
from pathlib import Path

def extract_hkl_from_folder():
    """Extract HKL parameter from current folder name"""
    current_dir = Path.cwd().name
    print(f"Current folder: {current_dir}")
    
    # Extract hkl from folder name (assuming format: HAP_XXX_template)
    parts = current_dir.split('_')
    if len(parts) >= 2:
        hkl_param = parts[1]
        if hkl_param.isdigit() and len(hkl_param) == 3:
            return hkl_param
    
    raise ValueError(f"Invalid folder name format. Expected HAP_XXX_template, got {current_dir}")

def generate_cp2k_input(hkl_param):
    """Generate CP2K input file with dynamic HKL parameter"""
    template_file = "cp2k_run/Template-HAP_hkl_md.inp"
    output_file = f"cp2k_run/HAP_{hkl_param}_md.inp"
    
    if os.path.exists(template_file):
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        content = content.replace("HAP_hkl", f"HAP_{hkl_param}")
        
        with open(output_file, 'w') as f:
            f.write(content)
        
        print(f"Generated CP2K input: {output_file}")
    else:
        print(f"Template file not found: {template_file}")

def generate_docker_compose(hkl_param):
    """Generate docker-compose.yml with dynamic HKL parameter"""
    template_file = "Template-docker-compose.yml"
    output_file = "docker-compose.yml"
    
    if os.path.exists(template_file):
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        content = content.replace("HAP_hkl", f"HAP_{hkl_param}")
        
        with open(output_file, 'w') as f:
            f.write(content)
        
        print(f"Generated docker-compose.yml")
    else:
        print(f"Template file not found: {template_file}")

def generate_lammps_input(hkl_param):
    """Generate LAMMPS input file with dynamic HKL parameter"""
    template_file = "lammps_run/Template-in.deepmd.lammps"
    output_file = "lammps_run/in.deepmd.lammps"
    
    if os.path.exists(template_file):
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Check if CO2 version exists
        co2_file = f"lammps_run/hap_{hkl_param}_co2.data"
        if os.path.exists(co2_file):
            # Use CO2 version
            content = content.replace("hap_hkl.data", f"hap_{hkl_param}_co2.data")
            print(f"Using CO2 version: hap_{hkl_param}_co2.data")
        else:
            # Use regular version
            content = content.replace("hap_hkl.data", f"hap_{hkl_param}.data")
            print(f"Using regular version: hap_{hkl_param}.data")
        
        with open(output_file, 'w') as f:
            f.write(content)
        
        print(f"Generated LAMMPS input: {output_file}")
    else:
        print(f"Template file not found: {template_file}")

def generate_specification_file(hkl_param):
    """Generate specification file with dynamic HKL parameter"""
    template_file = "Template-HAP_hkl_specification.txt"
    output_file = f"HAP_{hkl_param}_specification.txt"
    
    if os.path.exists(template_file):
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        content = content.replace("HAP_hkl", f"HAP_{hkl_param}")
        content = content.replace("hkl", hkl_param)
        
        with open(output_file, 'w') as f:
            f.write(content)
        
        print(f"Generated specification file: {output_file}")
    else:
        print(f"Template file not found: {template_file}")

def main():
    try:
        hkl_param = extract_hkl_from_folder()
        print(f"Extracted HKL parameter: {hkl_param}")
        
        # Generate configuration files
        generate_cp2k_input(hkl_param)
        generate_docker_compose(hkl_param)
        generate_lammps_input(hkl_param)
        generate_specification_file(hkl_param)
        
        print("Dynamic configuration generation completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
