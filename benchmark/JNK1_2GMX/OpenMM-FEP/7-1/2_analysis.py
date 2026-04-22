#!/usr/bin/env python

import sys
import os
import argparse
import numpy as np

from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer

sys.path.insert(0, "openmm_scripts")

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--work-dir', default='.',
                    help='Path to the working directory (default: current directory)')
parser.add_argument('-s', '--system', nargs=2, type=str, default=['complex', 'ligand'],
                    help='AFES system list (default: LIST["complex", "ligand"])')
parser.add_argument('-t', '--temperature', type=float, default=300.0,
                    help='Temperature (default: 300.0 K)')
parser.add_argument('-o', '--outfile', default='results.txt',
                    help='Output file name to store the calculated ΔΔG (default: results.txt)')
args = parser.parse_args()

T = args.temperature
kT_kcal = 0.0019872041 * T

dG = list()
outfile = open(args.outfile, 'w')
for system in args.system:
    # ========== 1. Preparation ==========
    target_dir = os.path.join(args.work_dir, system)
    nc_file = os.path.join(target_dir, f'openmm_fep_{system}.nc')
    reporter = MultiStateReporter(storage=nc_file, open_mode="r")

    # ========== 2. Calculate dG ==========
    analyzer = ReplicaExchangeAnalyzer(reporter)
    delta_f, d_delta_f = analyzer.get_free_energy()
    DeltaG_kcal = delta_f * kT_kcal
    dDeltaG_kcal = d_delta_f * kT_kcal
    dG.append(DeltaG_kcal[0, -1])
    outfile.write(f"{system} ΔG = {DeltaG_kcal[0, -1]:.4f} ± {dDeltaG_kcal[0, -1]:.4f} kcal/mol\n")

    # ========== 3. Per-lambda state overlap diagnostics ==========
    overlap = analyzer.mbar.compute_overlap()
    n_states = overlap['matrix'].shape[0]
    overlap_path = os.path.join(target_dir, f'openmm_fep_{system}_overlap_matrix.dat')
    with open(overlap_path, 'w') as overlap_file:
        overlap_file.write("# i   j   overlap\n")
        for i in range(n_states):
            for j in range(n_states):
                overlap_file.write(f"{i:3d} {j:3d} {overlap['matrix'][i, j]:12.6f}\n")

ddG = dG[0] - dG[1]
outfile.write(f"ΔΔG = {ddG:.4f} kcal/mol\n")
outfile.close()

