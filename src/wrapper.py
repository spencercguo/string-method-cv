"""
Wrapper script for string method calculations
"""
import numpy as np
import os
import os.system.chdir as cd
from subprocess import call
import sys
import argparse
import time
import asyncio
import pickle as pckl
from collections import namedtuple

from fileutils import *
from string import String

Params = namedtuple(
    'Params', 
    [
        'dt',           # Delta t for forward Euler calculation
        'smooth',       # smoothing parameter (s) 
        'k',            # spring constant for biased simulations (nonequilibration)
        't_us',         # umbrella sampling time used to calculate force vector
        't_us_eq',      # time to equilibrate umbrella 
        't_us_throw',   # time used in US that is not used for calculations
        't_eq',         # time for dragging umbrella
        't_eq_r',       # running equilibration time with umbrella
        'n_iter',       # number of string iterations to run 
        'timestep',     # MD timestep in units of 
        'stride',       # PLUMED stride 
        'tol',          # tolerance for iterations of string
        'keq'           # spring constant for equilibration
    ]
)

def parse():
    # Parse arguments to script
    parser = argparse.ArgumentParser(description="Run the string method in GROMACS")
    parser.add_argument(
        "--n", 
        type=int, 
        default=28, 
        dest=ncores, 
        help="Number of cores to run with."
    )
    parser.add_argument(
        "--name",
        type=str,
        dest=name,
        help="Name of the job and output files."
    )
    parser.add_argument(
        "-i",
        type=str,
        dest=initframes,
        help='Initial values of CVs for each image.'
    )
    parser.add_argument(
        "--cvfile",
        type=str,
        dest=cvfile,
        help="Name of file with CVs (.npy) of shape (R images, N CVs)."
    )
    parser.add_argument(
        "--folder",
        type=str,
        dest=folder,
        help="Name of folder header to write output files."
    )
    parser.add_argument(
        "--parfile",
        type=str,
        dest=parfile,
        help="Name of file with string run parameters."
    )
    nodetype=sys.argv[7]

    return ncores, names, initframes, cvfile, folder, parfile

def load_parfile(parfile):
    # loads parameter file which has inputs on one line
    # separated by commas
    with open(parfile, mode='r+') as f:
        params = f.readline().split(',')
    p = Param._make(params)

    return p

def umbrella_sample(i, z_target, cv_info, par):
    """
    Run dragging, equilibration, and restrained simulation on replicas of a 
    given string iteration.

    Parameters
    ---------
    i (int): the number of the iteration
    z_target (np.ndarray(R, n_CVs)): target CV values 
    """
    k = np.ones(len(z_curr)) * par.k
    keq = np.ones(len(z_curr)) * par.keq

    # compute or load current CV values

    plumed_outputs = [] 
    for m in range(1, NREPLICAS - 1):
        os.mkdir(f"im_{str(m).zfill(zim)}")
        cd(f"im_{str(m).zfill(zim)}")
        plumed_outputs.append("im_"+str(m).zfill(zim)+"/cv.txt")
        plumed_equil_file = fileutils.write_plumed_equil(
            cv_info,
            z_curr,
            z_target,
            par.t_eq,
            par.t_eq_r,
            par.k,
            par.keq,
            timestep=par.timestep,
            stride=10,
            whole=""
            output='plumed_equil.dat'
        )
        sb_equil = fileutils.write_sb_equil(
            cv_info,

        )
        sb_files.append(sb_equil)
        notphys.append(m)
        cd("..")

    # DO A BUNCH OF ASYNCIO STUFF TO RUN ALL THE IMAGES equilibration
    asyncio.run(run_all_images(sb_files)

    # check all images to ensure that within tolerance
    while notphysph:
        notphysph=[]
        for m in notphys:
            avg_z = np.mean(np.loadtxt(plumed_outputs[m]))[:,1:],axis=0)
            devs = [norm(cv-z) for z in zcurrent]
            close=np.argmin(devs)
            if close!=m:
                diff=np.abs(zcurrent[m]-cv)-np.abs(zcurrent[close]-cv)
                notphysph.append(m)
                kadd=np.asarray([int(x>-par.tol)*keq for x in diff])
                Keqnew[m]=Keq[m]+kadd
        # rerun equilibration if necessary
        for m in notphysph:
            cd(f"im_{str(m).zfill(zim)}")
            plumed_outputs.append("im_"+str(m).zfill(zim)+"/cv.txt")
            plumed_equil_file = fileutils.write_plumed_equil(
                cv_info,
                z_curr,
                z_target,
                par.t_eq,
                par.t_eq_r,
                par.k,
                par.keq,
                timestep=par.timestep,
                stride=10,
                whole=""
                output='plumed_equil.dat'
            )
            sb_equil = fileutils.write_sb_equil(
                cv_info,

            )
            sb_files.append(sb_equil)
            notphys.append(m)
            cd("..")
    
 
    # update z_curr
    z_curr = z_target
    # perform stationary umbrellas
    cv_output_files = []
    for m in range(1, NREPLICAS - 1):
        cd(f"im_{str(m).zfill(zim)}")
        plumed_umbr_file = fileutils.write_plumed_umbr(
            cv_info,
            z_curr,
            par.t_us_eq,
            par.t_us_throw,
            par.t_us,
            k,
            keq,
            timestep=par.timestep,
            stride=10,
            whole=""
        )
        sb_equil = fileutils.write_sb_equil(
            cv_info,

        )
        cv_output_files.append()
        notphys.append(m)
        cd("..")

    # DO A BUNCH OF STUFF TO LOAD CVS FROM PLUMED FILES
    z_curr = np.empty(
    for f in cv_output_files:


    return z_curr, derivs

async def run_all_images(sb_files):
    tasks = []
    for sb in sb_files:
        tasks.append(run_sb(sb))
    await asyncio.gather(tasks)

async def run_sb(sb_file):
    await call(f"sbatch {sb_file}")

def main():
    ncores, name, initframes, cvfile, folder, parfile = parse()

    print(f"""Settings\n
            ================================
            Number of cores:    {ncores}\n
            Name of job:        {name}\n
            Initial frame folder:{initframes}\n
            CV file:            {cvfile}\n
            Folder:             {folder}\n
            Parameter file:     {parfile}\n""")
    print("=================================\n")

    # Initialization steps
    print(f"Moving to {os.path.abspath(folder)}...")
    cd(folder)
    print(f"Making folder {os.path.abspath(name)}...")
    if not os.path.exists(name):
        os.mkdir(name) #edit this so that if name exists, it makes name_1, and so on.
    else:
        os.mkdir(f"{name}_1")
        name = f"{name}_1"
    print(f"Moving to {os.path.abspath(name)}...\n")
    cd(name)
    print("=================================\n")
    
    print(f"Initializing string...")
    # initialize string
    cv_info = infer_cvs(cvfile)
    z_target= np.load(cvfile)
    z_curr = np.zeros_like(z_curr)
    s = String(
        z_curr, dt=par.dt, k=par.k, smooth=par.smooth)
 
    par = load_parfile(parfile)

    n_images = s.R
    
    # padding zeros
    zeropad_im = 1 + int(np.log10(s.R))
    zeropad_it = 1 + int(np.log10(par.num_iter))
    
    # make data directories and copy initial frames
    os.mkdir(f"string_{str(0).zfill(zeropad_it)}")
    call(f"cp ../{initframes}/* string_{str(0).zfill(zit)}/", shell=True)

    #creating mdp files
    call("cp ../../nvt.mdp nvt_eq.mdp", shell=True)
    call("cp ../../nvt.mdp nvt_umbr.mdp", shell=True)
    
    with open("nvt_eq.mdp","a") as f:
        f.write("nsteps = "+str(int((tequil+tequilr)*1000/timestep))+"\n")
        f.write("dt = "+str(timestep))
        f.close()
    with open("nvt_umbr.mdp","a") as f:
        f.write("nsteps = "+str(int((tumbr+tumbr_throw+umbr_tequil)*1000/timestep))+"\n")
        f.write("dt = "+str(timestep))
        f.close()

    sf.makeplumedinit(CVinfo,num_CV)

    # z_target are the newly propagated CV values from the PREVIOUS ITERATION. 
    # Images need to be equilibrated to this point before you can do anything else.
    # z_curr are placeholders for values of CV just before zcurrent.
    # z_phys should be the physical values of CV after any set of calculations.
    # The above are currently updated only at the end of an iteration.

    #This is where the string iteration begins.

    print(f"Starting string iteration at {time.localtime(time.time())}...")
    k = np.ones_like(z_curr) * par.k
    keq = np.ones_like(z_curr) * par.keq

    for i in range(int(par.n_iter)):
        start = time.time()
        print(f"Iteration {i}")
        iter_folder = f"string_{i.zfill(zeropad_it)}"
        os.mkdir(iter_folder)
        cd(iter_folder)
        # perform dragging and equilibration for US
        cvs, deriv = umbrella_sample(i, z_curr)
        # Update string from simulation data
        s.update(cvs, deriv, par.smooth)
        s.write_to_file(iter_folder)
        cd("..")
        end = time.time()
        print(f"Iteration {i} took {end - start} seconds")
        print(f"======================================\n")
    
    print(f"String iteration complete!")

if __name__ == "__main__":
    main()
