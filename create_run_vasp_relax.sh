#!/bin/bash



cat > run_vasp_relax.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)



#################################
######     Relax       ##########
#################################

# ENCUTの値をもとに、INCARファイルを作成する
ENCUT=\$(tail -n -1 encut_cutoff.dat | awk '{print \$1}')

echo "relax structure INCAR" > INCAR
echo "ISTART = 0" >> INCAR
echo "ICHARG = 1" >> INCAR
echo "ISPIN = 1" >> INCAR
echo "ENCUT = \$ENCUT" >> INCAR
echo "EDIFF = 1.0e-4" >> INCAR
echo "LWAVE = .F." >> INCAR
echo "ISMEAR = -5" >> INCAR
echo "SIGMA = 0.01" >> INCAR
echo "PREC = accurate" >> INCAR
echo "LREAL = .F." >> INCAR
echo "NELM = 200" >> INCAR
echo "IBRION = 2" >> INCAR
echo "NSW = 1000" >> INCAR
echo "EDIFFG = -1.0e-2" >> INCAR
echo "ISIF = 3" >> INCAR
echo "NPAR = 4" >> INCAR

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4


hostname
date

EOF
echo "run_vasp_relax.sh file has been successfully created"

