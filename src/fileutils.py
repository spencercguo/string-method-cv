"""
Utility functions for file writing and reading
"""
import numpy as np
import os
import sys
from subprocess import call

def makeplumedinit(CVinfo,num_CV):
# writing first plumed file just to get first set of CVs

    f=open("plumedfileinit.dat","w+")
    for it in range(num_CV):
        CVdict[CVinfo[it][1]](CVinfo[it],f)

    f.write("PRINT STRIDE=1 ARG=* FILE=cv.txt")
    f.close()

def infer_cvs(cvfile):
    """
    """
    #File format name, type, array of atoms involved
    cv_info=[]
    with open(cvfile, mode='r+') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        cv_info.append(line.split()[:-1])
        cv_info[i].append([x for x in np.fromstring(line.split()[-1],sep=",")])

    return cv_info



def CVout(CVinfo,num_CV,num_images,zim,zit):
# writing first plumed file just to get first set of CVs

    zphys=np.zeros((num_images,num_CV))

    for m in range(1,num_images-1):
        oss("mpirun -np 1 plumed driver --plumed ../plumedfileinit.dat --igro "+str(m).zfill(zim)+".gro &>> ../out.txt ") #fix to zfill later to fit convention
        zphys[m,:]=np.genfromtxt("cv.txt")[1:]
        oss("rm cv.txt")

    return zphys

def write_plumed_equil(
    CVinfo,
    z_curr,
    z_target,
    t_eq,
    t_eq_r,
    k,
    keq,
    timestep=2, # FIX
    whole=None,
    stride=10,
    outfile='plumed.dat'
):
    """
    Writes an equilibration file in PLUMED to first drag CVs to their target
    values and then equilibration for a set amount of time.

    Parameters
    ---------
        z_curr (np.ndarray(N)): array of current CV values
        z_target (np.ndarray(N)): array of target CV values to drag to
        t_eq (int): how much time in # FIXXX to perform initial dragging
        t_eq_r (int): how much time in # to run umbrella equilibration
        k (int): force constant for initial dragging
        keq (int): force constant for umbrella equilibration
        timestep (int): MD timestep
        whole (str or None): for the WHOLEMOLECULES flag in PLUMED
        stride (int): PLUMED write stride
        outfile (str): name of the output PLUMED file
    """
    with open(outfile, mode='w+') as f:
        if whole is not None:
            f.write(f"WHOLEMOLECULES ENTITY0={whole} \n")

        for it in range(len(CVinfo)):
            CVdict[CVinfo[it][1]](CVinfo[it],f)
            f.write("restraint-%s:  ...\n" % CVinfo[it][0])
            f.write("    MOVINGRESTRAINT\n")
            f.write("    ARG=%s\n" % (CVinfo[it][0]))
            f.write("    AT0=%f STEP0=0 KAPPA0=%i \n" % (z_curr[it],k[it]))
            f.write("    AT1=%f STEP1=%i KAPPA1=%i \n" % (z_target[it],t_eq*1000/timestep,keq[it]))
            f.write("    AT2=%f STEP2=%i KAPPA2=%i \n" % (z_target[it], (t_eq + t_eq_r)*1000/timestep,keq[it]))
            f.write("...\n")

        f.write(f"PRINT STRIDE={stride} ARG=")
        for it in range(len(CVinfo)-1):
            f.write("%s," % CVinfo[it][0])
        f.write("%s FILE=cv.txt UPDATE_FROM=%f\n" % (CVinfo[len(CVinfo)-1][0],float(tequil1+tequil2/2)*1000))#*1000/timestep))
    return outfile

def write_plumed_umbr(
    CVinfo,
    z_curr,
    t_us_eq,
    t_us_throw,
    t_us,
    k,
    keq,
    timestep=0.002, # ps
    stride=10,
    whole=0
    outfile='plumed.dat'
):
    with open(outfile, mode='w+') as f:
        if whole!=0:
            f.write("WHOLEMOLECULES ENTITY0=%s \n" % (whole))

        for it in range(len(CVinfo)):
            CVdict[CVinfo[it][1]](CVinfo[it],f)
            f.write(f"restraint-%s:  ...\n" % CVinfo[it][0])
            f.write(f"    MOVINGRESTRAINT\n")
            f.write(f"    ARG=%s\n" % (CVinfo[it][0]))
            f.write(f"    AT0={z_curr[it]} STEP0=0 KAPPA0={keq[it]} \n")
            step1 = t_us_eq * 1000 / timestep
            f.write(f"    AT1={z_curr[it]} STEP1={step1} KAPPA1={k[it]} \n")
            step2 = (t_us_eq + t_us_throw + t_us)*1000/timestep,
            f.write(f"    AT2={z_curr[it]} STEP2={step2} KAPPA2={k[it]}\n")
            f.write("...\n")

        f.write(f"PRINT STRIDE={stride} ARG=")
        for it in range(len(CVinfo)-1):
            f.write("%s," % CVinfo[it][0])
        f.write("%s FILE=cv.txt UPDATE_FROM=%i\n" % (CVinfo[len(CVinfo)-1][0],(umbr_tequil+tumbr_throw)*1000))#*/timestep))

        for it in range(len(CVinfo)):
            f.write("DUMPDERIVATIVES STRIDE=%i ARG=%s FILE=deriv_%i UPDATE_FROM=%f\n" % (STRIDE,CVinfo[it][0],it,(umbr_tequil+tumbr_throw)*1000))#*1000/timestep))



def write_gromacs_sb(
    i,
    m,
    mdp,
    groc,
    topp,
    deffnm,
    cnp,
    plmd,
    time,
    zim=2,
    FwdE=0,
    nodetype="cc",
    outfile="run.sbatch"
):
    with open(outfile, mode='w+') as f:
        name = f"i_{i}_im_{m}_{deffnm}"
        if groc == "equil" and deffnm == "equil":
            if nodetype[-2:]=="sb":
                cnp=cnp*4
            else:
                if cnp==1:
                    cnp=4
                else:
                    cnp=28
        f.write("""
            #!/bin/bash 
            #SBATCH --job-name={name}
            #SBATCH --output={name}.out
            #SBATCH --error={name}.err""")
        if nodetype=="cc":
            f.write("""#SBATCH --partition=broadwl
                    #SBATCH --account=pi-dinner""")
        elif nodetype=="wd":
            f.write("#SBATCH --partition=weare-dinner2 \n#SBATCH --qos=weare-dinner --account=weare-dinner \n")
        elif nodetype=="ccsb":
            f.write("#SBATCH --partition=sandyb \n#SBATCH --account=pi-dinner \n")
        elif nodetype=="wdsb":
            f.write("#SBATCH --partition=weare-dinner1 \n#SBATCH --qos=weare-dinner --account=weare-dinner \n")

        if cnp%28==0:
            f.write("#SBATCH --ntasks=%i\n#SBATCH --exclusive\n#SBATCH --time=%i:00 \n" % (cnp,min(int(time*200/np.sqrt(cnp))+1,2160) ))
        elif cnp%16==0:
            f.write("#SBATCH --ntasks=%i\n#SBATCH --exclusive\n#SBATCH --time=%i:00 \n" % (cnp,min(int(time*200/np.sqrt(cnp))+1,2160) ))
        else:
            f.write("#SBATCH --ntasks=%i --nodes=1\n#SBATCH --time=%i:00 \n" % (cnp,min(int(time*150/np.sqrt(cnp))+1,2160) )) # --cpus-per-task=%i



        if nodetype[-2:]=="sb":
            f.write("module unload intel\nmodule unload intelmpi\nmodule unload openmpi \nmodule load intelmpi/5.0+intel-15.0 \nmodule load python/2.7-2015q2\n")
        else:
            f.write("module purge \nmodule load gromacs \nmodule load plumed \nmodule load python\n")
        f.write("mpirun -np 1 gmx_mpi grompp -f %s.mdp -c %s.gro -p %s.top -o %s.tpr -maxwarn 2\n"\
                %(mdpf,groc,topp,deffnm))
        f.write("mpirun -np %i gmx_mpi mdrun -deffnm %s -plumed %s.dat -ntomp 1\n"\
                %(cnp,deffnm,plmd))
        if FwdE==1:
            f.write("python ../../../../stringstep1.py im_"+str(m).zfill(zim)+" "+str(m))


