#!/usr/bin/env python3

import os
import subprocess

# Save original PATH and LD_LIBRARY_PATH
original_path = os.environ.get('PATH', '')
original_ld_library_path = os.environ.get('LD_LIBRARY_PATH', '')

# Source environment.sh and capture the modified environment
# Use a shell command to source environment.sh and print PATH and LD_LIBRARY_PATH
command = f"source $ISSM_DIR/etc/environment.sh && echo $PATH && echo $LD_LIBRARY_PATH"
result = subprocess.run(command, shell=True, executable="/bin/bash", capture_output=True, text=True, check=True)

# Split the output into PATH and LD_LIBRARY_PATH
output_lines = result.stdout.strip().split('\n')
if len(output_lines) >= 2:
    new_path = output_lines[0]
    new_ld_library_path = output_lines[1]
else:
    raise RuntimeError("Failed to capture PATH and LD_LIBRARY_PATH from environment.sh")

# Prepend original paths to ensure mvapich2 and gcc/12 take precedence
os.environ['PATH'] = f"{original_path}:{new_path}"
os.environ['LD_LIBRARY_PATH'] = f"{original_ld_library_path}:{new_ld_library_path}"

# Run MATLAB with the input file
subprocess.run(["matlab", "-nodesktop", "-nosplash", "-r", "run('issm_env'); run('runme.m'); exit"], check=True)