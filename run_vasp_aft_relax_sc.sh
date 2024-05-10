#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd $PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=$(cat $PBS_NODEFILE | wc -l)






############################################
###   final SC calculation after relax   ###
############################################


# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p final_sc
mkdir -p bandgap_cal

cp POTCAR  ./final_sc
    cp POSCAR_gt ./final_sc/POSCAR
    cp POSCAR_gt ./bandgap_cal/POSCAR

cd final_sc

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_ENCUT=$(cat ../BEST_ENCUT.dat)
echo "final sc-calculation INCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "$BEST_ENCUT" >> INCAR
echo "ISMEAR = -5; SIGMA = 0.1" >> INCAR
echo "PREC = accurate" >> INCAR
echo "ISPIN = 2" >> INCAR
echo "EDIFF = 1.0e-6 # default:10^-4, SCF計算の収束条件,推奨は1E-6らしい。" >> INCAR
cat ../INCAR_tail >> INCAR


# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_KPOINTS=$(cat ../BEST_KPOINTS.dat)
echo "BEST KPOINTS"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "$BEST_KPOINTS" >> KPOINTS
echo "0 0 0" >> KPOINTS





echo "################################################"
echo "################################################"
echo "###  VASP calculation for final SC starts!   ###"
echo "################################################"
echo "################################################"
echo ""

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
date
echo "#######################################"
echo "###            VASP ends!           ###"
echo "#######################################"


cp CHGCAR ../bandgap_cal/

cd ../

