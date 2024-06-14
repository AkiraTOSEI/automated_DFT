#!/bin/bash


# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <RELAX>"
    exit 1
fi
RELAX="$1"

cat > run_vasp_relax.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)

EOF







if [ "$RELAX" = "NO" ]; then
    cat >> run_vasp_relax.sh <<EOF
echo " NO relaxation !!!!!!!! "
EOF
    echo "run_vasp_relax.sh file has been successfully created"

elif [ "$RELAX" = "YES" ]; then
    cat >> run_vasp_relax.sh <<EOF
#################################
######     Relax       ##########
#################################

mkdir -p relax_cal

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
python calculate_EDIFFG.py > EDIFFG.dat
python calculate_EDIFF.py > EDIFF.dat
cp POSCAR_gt ./relax_cal/POSCAR
cp POTCAR ./relax_cal
cd ./relax_cal

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_KPOINTS=\$(cat ../BEST_KPOINTS.dat)
MESH_METHOD=\$(cat ../mesh_method.dat)
echo "BEST KPOINTS"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "\$MESH_METHOD" >> KPOINTS
echo "\$BEST_KPOINTS" >> KPOINTS
echo "0 0 0" >> KPOINTS

# ENCUTの値をもとに、INCARファイルを作成する
BEST_ENCUT=\$(cat ../BEST_ENCUT.dat)
EDIFF=\$(cat ../EDIFF.dat)
EDIFGF=\$(cat ../EDIFFG.dat)
echo "relax structure INCAR" > INCAR
echo "ISTART = 0" >> INCAR
echo "ICHARG = 1" >> INCAR
echo "ISPIN = 2" >> INCAR
echo "\$BEST_ENCUT" >> INCAR
echo "\$EDIFF" >> INCAR
echo "LWAVE = .F." >> INCAR
echo "ISMEAR = 0" >> INCAR
echo "SIGMA = 0.03" >> INCAR
echo "PREC = accurate" >> INCAR
echo "LREAL = .F." >> INCAR
echo "NELM = 200" >> INCAR
echo "IBRION = 2" >> INCAR
echo "NSW = 1000" >> INCAR
echo "\$EDIFF" >> INCAR
echo "ISIF = 3" >> INCAR
echo "NPAR = 4" >> INCAR
echo "NCORE = 8 #" >> INCAR
echo "KPAR = 4 # defalut:コア数、並列計算に使うコア数" >> INCAR
cat ../incar_magmom.dat >> INCAR
cat ../INCAR_tail >> INCAR

echo "################################################"
echo "################################################"
echo "###  VASP calculation for  Relax starts!     ###"
echo "################################################"
echo "################################################"
echo ""

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4

cd ../

hostname
date

echo "#######################################"
echo "###           Relax Ends!           ###"
echo "#######################################"
EOF
    echo "run_vasp_relax.sh file has been successfully created"

else
    echo "<RELAX> argment error!!!!!!!! RELAX:$RELAX "
fi




