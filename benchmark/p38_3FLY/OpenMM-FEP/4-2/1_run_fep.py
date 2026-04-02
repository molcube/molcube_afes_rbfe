#!/usr/bin/env python

import os
import subprocess
import threading
import argparse
from datetime import datetime


def log_message(message):
    print(f"{message} - {datetime.now()}")
    return


class SLURM:
    def __init__(self, job_name, work_dir='.', target_dir='.', ncpu=1, ngpu=0, output_prefix='slurm-%j'):
        self.job_name = job_name
        self.work_dir = work_dir
        self.target_dir = target_dir
        self.ncpu = ncpu
        self.ngpu = ngpu
        self.output_prefix=output_prefix
        self.script_dir  = os.path.abspath(os.path.join(self.work_dir, 'job_scripts'))
        self.script_file = os.path.join(self.script_dir, f"{self.job_name}.sh")

    def get_header(self):
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --job-name={self.job_name}\n"
            f"#SBATCH --output={self.target_dir}/{self.output_prefix}.out\n"
            f"#SBATCH --error={self.target_dir}/{self.output_prefix}.err\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={self.ncpu}\n"
            f"#SBATCH --exclude=node[001-002]\n"
        )
        if self.ngpu > 0:
            header += (
                f"#SBATCH --gres=gpu:{self.ngpu}\n\n"
                f"CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS\n"
            )
        return header
    def write(self, cmd):
        os.makedirs(self.script_dir, exist_ok=True)
        header = self.get_header()

        start_cmd = """eval "$(/home/apps/miniconda3/bin/conda shell.bash hook)"
    conda activate openmm
    start=$(date +%s)
    """

        end_cmd = """end=$(date +%s)
    echo "Total time: $((end - start)) seconds"
    """

        with open(self.script_file, "w") as f:
            f.write(header)
            f.write("\n")
            f.write(start_cmd)
            f.write("\n")
            f.write(cmd)
            f.write("\n")
            f.write(end_cmd)
    def submit(self, dependency=None):
        cmd = ["sbatch"]
        if dependency:
            cmd += ["--dependency=afterok:" + dependency]
        cmd.append(self.script_file)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_id = result.stdout.strip().split()[-1]
        log_message(f"Submitted job {self.job_name} with ID {job_id}")
        return job_id


class BASH:
    def __init__(self, job_name, work_dir='.', target_dir='.', output_prefix='output'):
        self.job_name = job_name
        self.work_dir = work_dir
        self.target_dir = target_dir
        self.output_prefix=output_prefix
        self.script_dir  = os.path.abspath(os.path.join(self.work_dir, 'job_scripts'))
        self.script_file = os.path.join(self.script_dir, f"{self.job_name}.sh")

    def get_header(self):
        header = f"#!/bin/bash\n"
        return header

    def write(self, cmd):
        os.makedirs(self.script_dir, exist_ok=True)
        header = self.get_header()
        with open(self.script_file, "w") as f:
            f.write(header)
            f.write("\n")
            f.write(cmd)

    def submit(self):
        cmd = ["bash"]
        cmd.append(self.script_file)

        outfile = open(os.path.join(self.target_dir, f'{self.output_prefix}.log'), 'w')
        errfile = open(os.path.join(self.target_dir, f'{self.output_prefix}.err'), 'w')

        process = subprocess.Popen(cmd, stdout=outfile, stderr=errfile, text=True)
        log_message(f"Submitted job {self.job_name} with ID {process.pid}")
        return process


class OpenMM_FEP:
    def __init__(self, base_dir='.'):
        self.base_dir = base_dir

    def get_cmd(self, target_dir, system):
        inpfile = f"{target_dir}/step2_run_fep.inp"
        sysfile = f"{target_dir}/step1_input.xml"
        pdbfile = f"{target_dir}/step1_input.pdb"

        cmd = f"python -u openmm_scripts/openmm_run.py -i {inpfile} -f {sysfile} -c {pdbfile}" \
              f" -o {target_dir}/openmm_fep_{system}"
        return cmd

    def run_steps(self, system, script_method='slurm'):
        target_dir = os.path.join(self.base_dir, system)

        job_name = f'openmm_fep_{system}'
        if script_method == 'slurm':
            script = SLURM(job_name, work_dir=self.base_dir, target_dir=target_dir, ncpu=1, ngpu=1)
        elif script_method == 'bash':
            script = BASH(job_name, work_dir=self.base_dir, target_dir=target_dir)
        cmd = self.get_cmd(target_dir, system)
        script.write(cmd)
        script.submit()

    def run(self, script_method='slurm'):
        systems = ['complex', 'ligand']
        threads = []
        for system in systems:
            t = threading.Thread(target=self.run_steps, args=(system, script_method))
            threads.append(t)

        for t in threads:
            t.start()

        for t in threads:
            t.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--work-dir', default='.',
                        help='Path to the working directory (default: current directory)')
    parser.add_argument('-m', '--method',   default='slurm', choices=['slurm', 'bash'],
                        help='Job submission method (default: slurm)')
    args = parser.parse_args()

    openmm_fep = OpenMM_FEP(args.work_dir)
    openmm_fep.run(args.method)

