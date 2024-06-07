#!/bin/bash





# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <RELAX>"
    exit 1
fi
RELAX="$1"

cat > run_vasp_aft_relax_sc.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)







echo "############################################"
echo "###   final SC calculation after relax   ###"
echo "############################################"



# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p final_sc
mkdir -p bandgap_cal

cp POTCAR  ./final_sc
EOF

if [ "$RELAX" = "YES" ]; then
    cat >> run_vasp_aft_relax_sc.sh <<EOF
    cp ./relax_cal/CONTCAR ./final_sc/POSCAR
    cp ./relax_cal/CONTCAR ./bandgap_cal/POSCAR
EOF
elif [ "$RELAX" = "NO" ]; then
    cat >> run_vasp_aft_relax_sc.sh <<EOF
    cp POSCAR_gt ./final_sc/POSCAR
    cp POSCAR_gt ./bandgap_cal/POSCAR
EOF
else
    echo "<RELAX> argment error!!!!!!!! "
fi


cat >> run_vasp_aft_relax_sc.sh <<EOF

python calculate_EDIFF.py > EDIFF.dat
cd final_sc

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_ENCUT=\$(cat ../BEST_ENCUT.dat)
EDIFF=\$(cat ../EDIFF.dat)
echo "final sc-calculation INCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "\$BEST_ENCUT" >> INCAR
echo "ISMEAR = -5; SIGMA = 0.1" >> INCAR
echo "PREC = accurate" >> INCAR
echo "ISPIN = 2" >> INCAR
echo "\$EDIFF" >> INCAR
cat ../INCAR_tail >> INCAR
cat ../incar_magmom.dat >> INCAR

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_KPOINTS=\$(cat ../BEST_KPOINTS.dat)
MESH_METHOD=\$(cat ../mesh_method.dat)
echo "BEST KPOINTS"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "\$MESH_METHOD" >> KPOINTS
echo "\$BEST_KPOINTS" >> KPOINTS
echo "0 0 0" >> KPOINTS





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


E_fermi=\$(grep "E-fermi" OUTCAR | awk '{print \$3}')
if [ -z "\$E_fermi"]; then
    echo "Gamma" > ../mesh_method.dat
    MESH_METHOD=\$(cat ../mesh_method.dat)
    echo "BEST KPOINTS with Gamma"
    echo "k-points" > KPOINTS
    echo "0" >> KPOINTS
    echo "\$MESH_METHOD" >> KPOINTS
    echo "\$BEST_KPOINTS" >> KPOINTS
    echo "0 0 0" >> KPOINTS
    
    echo "Failure to get  Fermi Enegry. Then re-calculate again with KPOINTS-Gamma"  
    echo "###################################################"
    echo "###################################################"
    echo "###  VASP calculation for final SC re-starts!   ###"
    echo "###################################################"
    echo "###################################################"
    echo ""
    
    mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
    date
    echo "#######################################"
    echo "###            VASP ends!           ###"
    echo "#######################################"
fi



cp CHGCAR ../bandgap_cal/

cd ../

EOF

echo "run_vasp_aft_relax_sc.sh file has been successfully created"

