#!/bin/bash

# Check if a file name is provided as the first argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <FILE_NAME>"
  exit 1
fi

FILE_NAME="$1"



# コマンドライン引数から.vaspファイルのパスを取得

cat > run_vasp_test_sc.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N $FILE_NAME
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)




echo "############################################"
echo "###   test SC calculation after relax   ###"
echo "############################################"
echo "Monkhorst Packで収束するかを確認するためのテスト計算"


# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p test_sc

cp POTCAR  ./test_sc
cp POSCAR_gt ./test_sc/POSCAR
EOF


cat >> run_vasp_test_sc.sh <<EOF

python calculate_EDIFF.py > EDIFF.dat

# test計算用のENCUT、KPOINTSを設定する
echo "ENCUT = 520" > BEST_ENCUT.dat
python create_kpoints.py > BEST_KPOINTS.dat
echo "Monkhorst Pack" > mesh_method.dat



# K点の対称性によって、対称性操作を切るかどうかを決める
read kx ky kz < ./BEST_KPOINTS.dat
if ! [[ "\$kx" =~ ^[0-9]+$ && "\$ky" =~ ^[0-9]+$ && "\$kz" =~ ^[0-9]+$ ]]; then
  echo "Error: The file does not contain three integers."
  exit 1
fi
# Determine if all three numbers are equal
if [ "\$kx" -eq "\$ky" ] && [ "\$ky" -eq "\$kz" ]; then
  # If all numbers are equal, create an empty file
  > INCAR_tail
else
  # If numbers are not equal, write specified content to the file
  echo "ISYM = -1" > INCAR_tail
  echo "SYMPREC = 0.00000001" >> INCAR_tail
fi

# test_sc ディレクトリに入ってテスト計算を行う
cd test_sc

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_ENCUT=\$(cat ../BEST_ENCUT.dat)
EDIFF=\$(cat ../EDIFF.dat)
echo "final sc-calculation INCAR" > INCAR
echo "ISTART = 0 ; ICHARG = 2" >> INCAR
echo "\$BEST_ENCUT" >> INCAR
echo "ISMEAR = -5; SIGMA = 0.1" >> INCAR
#echo "PREC = accurate" >> INCAR # テスト計算なので、accurateを使わない
echo "ISPIN = 2" >> _INCAR
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
echo "###  VASP calculation for test SC starts!   ###"
echo "################################################"
echo "################################################"
echo ""

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
date


E_fermi=\$(grep "E-fermi" OUTCAR | awk '{print \$3}')
if [ -z "\$E_fermi"]; then
    echo "Gamma" > ../mesh_method.dat
    echo "KPOINS -> Monkhorst Pack -> Gamma"
fi
echo "#######################################"
echo "###       test VASP ends!           ###"
echo "#######################################"




cd ../

EOF

echo "run_vasp_test_sc.sh file has been successfully created"

