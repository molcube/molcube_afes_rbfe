#!/usr/bin/env python

import os
import argparse

from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer

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
    target_dir = os.path.join(args.work_dir, system)
    nc_file = os.path.join(target_dir, f'openmm_fep_{system}.nc')

    # ========== 1. Open existing nc file ==========
    reporter = MultiStateReporter(storage=nc_file, open_mode="r")

    # ========== 2. Run analysis ==========
    analyzer = ReplicaExchangeAnalyzer(reporter)
    delta_f, d_delta_f = analyzer.get_free_energy()

    DeltaG_kcal = delta_f * kT_kcal
    dDeltaG_kcal = d_delta_f * kT_kcal

    dG.append(DeltaG_kcal[0][-1])
    outfile.write(f"{system} ΔG = {DeltaG_kcal[0][-1]:.4f} ± {dDeltaG_kcal[0][-1]:.4f} kcal/mol\n")

ddG = dG[0] - dG[1]
outfile.write(f"ΔΔG = {ddG:.4f} kcal/mol\n")

