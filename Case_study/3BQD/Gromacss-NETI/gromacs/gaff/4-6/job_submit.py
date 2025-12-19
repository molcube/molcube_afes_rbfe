import os
import subprocess
import multiprocessing
import argparse
from datetime import datetime

# mdp filenames
# ti_l1.mdp
# ti_l0.mdp
# eq_nvt_l1.mdp
# eq_nvt_l0.mdp
# eq_npt_l1.mdp
# eq_npt_l0.mdp
# em_l1.mdp
# em_l0.mdp

def log_message(message):
    print(f"{message} - {datetime.now()}")
    return

def submit_job(script_content, job_name, dependency=None):
    script_dir = os.path.abspath('slurm_scripts')
    os.makedirs(script_dir, exist_ok=True)
    script_file = os.path.join(script_dir, f"{job_name}.sh")

    with open(script_file, "w") as f:
        f.write(script_content)

    cmd = ["sbatch"]
    if dependency:
        cmd += ["--dependency=afterok:" + dependency]
    cmd.append(script_file)

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = result.stdout.strip().split()[-1]
    log_message(f"Submitted job {job_name} with ID {job_id}")
    return job_id

def create_slurm_script(target, mdp_file, output_prefix, job_prefix, system, stage,
                        crd_file=None, ncpu=6, ngpu=1,
                        additional_grompp_args="",
                        additional_mdrun_args=""):
    work_dir = os.path.join(os.path.abspath(f"{target}"), f"{stage}")
    mdp_file_abs = os.path.abspath(mdp_file)
    topol = os.path.abspath(f"{system}/topol.top")
    index = os.path.abspath(f"{system}/index.ndx")
    if crd_file == None:
        crd_file = os.path.abspath(f"{system}/step1_input.gro")

    grompp_cmd = f"gmx grompp -f {mdp_file_abs} -p {topol} -c {crd_file} -n {index} -r {crd_file}" \
                 f" -o {work_dir}/{output_prefix}.tpr -po {work_dir}/{output_prefix}out.mdp" \
                 f" -maxwarn 9999 {additional_grompp_args}"

    os.makedirs(os.path.join(work_dir, stage), exist_ok=True)

    script_content = f"""#!/bin/bash
#SBATCH --job-name={job_prefix}{system}_{stage}_{output_prefix}
#SBATCH --output={work_dir}/{output_prefix}.out
#SBATCH --error={work_dir}/{output_prefix}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={ncpu}
"""

    if ngpu > 0:
        script_content += f"""#SBATCH --gres=gpu:{ngpu}

CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS
"""

    script_content += f"""
module load gromacs
{grompp_cmd}
gmx mdrun -nobackup -deffnm {work_dir}/{output_prefix}
"""
    return script_content


def run_steps(target, state_num, job_prefix, system):
    # MINIMIZATION
    output_prefix = "em"
    em_script = create_slurm_script(target, f"mdp_files/em_l{state_num}.mdp", output_prefix, job_prefix, system, "em") 
    em_job_id = submit_job(em_script, f"em_{target.replace('/', '_')}")

    # NVT EQUILIBRATION
    previous_gro = os.path.abspath(f"{target}/em/em.gro")
    mdp_file = f"mdp_files/eq_nvt_l{state_num}.mdp"
    output_prefix = "eq_nvt"
    eq_script = create_slurm_script(target, mdp_file, output_prefix, job_prefix, system, "eq", previous_gro)
    eq_job_id = submit_job(eq_script, f"eq_nvt_{target.replace('/', '_')}", dependency=em_job_id)

    # NPT EQUILIBRATION
    previous_gro = os.path.abspath(f"{target}/eq/eq_nvt.gro")
    mdp_file = f"mdp_files/eq_npt_l{state_num}.mdp"
    output_prefix = "eq"
    eq_script = create_slurm_script(target, mdp_file, output_prefix, job_prefix, system, "eq", previous_gro)
    eq_job_id = submit_job(eq_script, f"eq_{target.replace('/', '_')}", dependency=eq_job_id)

    # TRANSITIONS
    t0_file = os.path.abspath(f"{target}/transitions/t0.txt")
    os.makedirs(os.path.dirname(t0_file), exist_ok=True)
    with open(t0_file, "w") as f:
        f.write("0 \n")

    eq_trr = os.path.abspath(f"{target}/eq/eq.trr")
    eq_tpr = os.path.abspath(f"{target}/eq/eq.tpr")
    frame_gro = os.path.abspath(f"{target}/transitions/frame.gro")

    trjconv_cmd = f"gmx trjconv -f {eq_trr} -s {eq_tpr} -o {frame_gro} -ur compact -pbc mol -sep -b 100 < {t0_file}"

    trjconv_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_prefix}trjconv_{target.replace('/', '_')}
#SBATCH --output={os.path.abspath(f"{target}/transitions/trjconv_{target.replace('/', '_')}.out")}
#SBATCH --error={os.path.abspath(f"{target}/transitions/trjconv_{target.replace('/', '_')}.err")}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gromacs
{trjconv_cmd}
rm {t0_file}
"""
    trjconv_job_id = submit_job(trjconv_script_content, f"trjconv_{target.replace('/', '_')}", dependency=eq_job_id)

    for i in range(0, 100, 10):
        job_indices = list(range(i, min(i + 10, 100)))
        output_prefix = f"ti_{i // 10}"
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_prefix}{system}_{state_num}_{output_prefix}
#SBATCH --output={os.path.abspath(f'{target}/transitions/{output_prefix}.out')}
#SBATCH --error={os.path.abspath(f'{target}/transitions/{output_prefix}.err')}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1

CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS

module load gromacs
"""
        for idx in job_indices:
            mdp_file = f"mdp_files/ti_l{state_num}.mdp"
            script_content += f"""
echo {idx}
date
gmx grompp -f {mdp_file} -n ./{system}/index.ndx -c ./{target}/transitions/frame{idx}.gro -r ./{target}/transitions/frame{idx}.gro -p ./{system}/topol.top -o ./{target}/transitions/ti{idx}.tpr -po ./{target}/transitions/ti{idx}.mdp -maxwarn 9999
gmx mdrun -s ./{target}/transitions/ti{idx}.tpr -e ./{target}/transitions/ti{idx}.edr -c ./{target}/transitions/ti{idx}.gro -o ./{target}/transitions/ti{idx}.trr -g ./{target}/transitions/ti{idx}.log -dhdl ./{target}/transitions/dhdl{idx}.xvg -x ./{target}/transitions/ti{idx}.xtc -cpo ./{target}/transitions/ti{idx}.cpt
"""
        submit_job(script_content, f"ti_{i // 10}_{target.replace('/', '_')}", dependency=trjconv_job_id)


def process_combination(job_prefix, system, state, trial):
    target = f"{system}/{state}/{trial}"
    os.makedirs(f"{target}/em", exist_ok=True)
    os.makedirs(f"{target}/eq", exist_ok=True)
    os.makedirs(f"{target}/transitions", exist_ok=True)
    state_num = 0 if state == "stateA" else 1

    run_steps(target, state_num, job_prefix, system)


def run_process_combination(combo):
    process_combination(*combo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--systems', type=str, nargs=2, default=['complex', 'ligand'],
                        help='Name of FEP Systems (default: ["complex", "ligand"]).')
    parser.add_argument('-r', '--replica', type=int, required=True, help='Number of Replica(s)')
    parser.add_argument('--job-prefix', type=str, default='', help='SLURM Job Name Prefix')
    args = parser.parse_args()

    log_message("start")

    job_prefix = args.job_prefix
    if job_prefix: job_prefix = job_prefix + '_'
    systems = args.systems
    states = ["stateA", "stateB"]
    trials = [f"run{i+1}" for i in range(args.replica)]

    combinations = [(job_prefix, system, state, trial) for system in systems for state in states for trial in trials]

    with multiprocessing.Pool(processes=len(combinations)) as pool:
        pool.map(run_process_combination, combinations)

    log_message("finish")

