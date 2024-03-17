#!/bin/bash






cat > run_vasp_aft_relax_sc.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)



############################################
###   final SC calculation after relax   ###
############################################


# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p final_sc
mkdir -p bandgap_cal
cp POTCAR KPOINTS ./final_sc
cp CONTCAR ./final_sc/POSCAR
cp CONTCAR ./bandgap_cal/POSCAR

cd final_sc

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
ENCUT=\$(tail -n -1 ../encut_cutoff.dat | awk '{print \$1}')
echo "final sc-calculation INCAR" > INCAR
ecoo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "ENCUT = \$ENCUT" >> INCAR
echo "ISMEAR = -5; SIGMA = 0.1" >> INCAR
echo "PREC = accurate" >> INCAR


echo "################################################"
echo "################################################"
echo "###  VASP calculation for final SC starts!   ###"
echo "################################################"
echo "################################################"
echo ""

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
date
echo "#######################################"
echo "###            VASP ends!           ###"
echo "#######################################"


cp CHGCAR ../bandgap_cal/

cd ../

EOF

echo "run_vasp_aft_relax_sc.sh file has been successfully created"

