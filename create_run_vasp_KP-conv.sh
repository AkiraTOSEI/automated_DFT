#!/bin/bash

# Check if a file name is provided as the first argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <FILE_NAME>"
  exit 1
fi

FILE_NAME="$1"

# Ensure kpoints_candidate.csv exists
if [ ! -f "kpoints_candidate.csv" ]; then
  echo "K-points candidates file was not generated."
  exit 1
fi

# Start creating run_vasp_KP-conv.sh
cat > run_vasp_KP-conv.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N $FILE_NAME
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)




echo "########################################################"
echo "########################################################"
echo "###  VASP calculation for KPOINTS convergence stars! ###"
echo "########################################################"
echo "########################################################"

cp POSCAR_distorted POSCAR

# 収束判定を初期化

EOF

i=1
while IFS=, read -r kx ky kz; do

    cat >> run_vasp_KP-conv.sh <<EOF
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: $kx $ky $kz
echo "Processing k-point $kx $ky $kz"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "$kx $ky $kz" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=\$(grep 'free  ' OUTCAR | tail -1 | awk '{print \$5}')
FORCE=\$(awk '/TOTAL-FORCE/ {getline; getline; print \$4, \$5, \$6}' OUTCAR) 
PRESSURE=\$(grep "external pressure" OUTCAR | awk '{print \$4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "\${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "\${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "\${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

EOF

    if [ $i -eq 1 ]; then
        cat >> run_vasp_KP-conv.sh <<EOF
echo "$kx $ky $kz"  > kpoints.dat
echo "\$energy" > energy_KP.dat
echo \$FORCE > force_KP.dat 
echo \$PRESSURE > pressure_KP.dat
EOF
    else
        cat >> run_vasp_KP-conv.sh <<EOF
echo "$kx $ky $kz" >> kpoints.dat
echo "\$energy" >> energy_KP.dat
echo \$FORCE >> force_KP.dat 
echo \$PRESSURE >> pressure_KP.dat
EOF
    fi
    i=$((i+1))

done < kpoints_candidate.csv

cat >> run_vasp_KP-conv.sh <<EOF


#収束判定ができるcppファイルをコンパイルする 
g++ -std=c++11 -o check_conv convergence_checker.cpp
g++ -std=c++11 -o calculate_l2_norm l2_norm.cpp
g++ -std=c++11 -o calculate_vol calculate_vol.cpp
g++ -std=c++11 -o convergence_summary convergence_summary.cpp

# 原子数を取得する
total_lines=\$(wc -l < POSCAR)
NUM_ATOM=\$((total_lines - 8))

# Pressureの最終値との差分が5 GPa, 絶対値が5 GPa以下になるようにする. 
# VASPでは圧力単位がKilo Barなので、GPaにへんかんするために0.1をかける
./check_conv pressure_KP.dat 0.1 0.5 last > pressure_diff.dat
./check_conv pressure_KP.dat 0.1 5. zero > pressure_abs.dat
# エネルギー/原子の最終値との差分が1meV以下になるようにする
./check_conv energy_KP.dat "1000/\$NUM_ATOM" 1 last > energy_KP_diff.dat
# 力はl2_normをとり、その値が
./calculate_l2_norm force_KP.dat > force_norm.dat
./check_conv force_norm.dat 1. 0.01 last > force_diff.dat

# KPOINTの体積を計算する.
./calculate_vol kpoints.dat > k_vol.dat

# 結果ファイルを作成する
paste force_diff.dat pressure_abs.dat pressure_diff.dat energy_KP_diff.dat > __tmp_conv.dat
echo "Fx Fy Fz F_norm(eV/Å) F_diff(eV/Å) pressure(GPa) P_diff(GPa) Energy(meV/atom) E_diff(meV/atom)" > __kp_history.dat
awk '{print \$1, \$2, \$7, \$8, \$10, \$11}' __tmp_conv.dat > tmp.dat; paste force_KP.dat tmp.dat| expand | tr -s "[:space:]" >> __kp_history.dat
echo  "kx ky kz k_vol" > kp_tmp1.dat
paste kpoints.dat k_vol.dat > kp_tmp2.dat 
cat kp_tmp1.dat kp_tmp2.dat > kp_tmp.dat
paste kp_tmp.dat __kp_history.dat | expand | tr -s "[:space:]" >  kp_history.dat
# 収束判定のファイルを作成する
awk '{print \$3, \$6, \$9, \$12}' __tmp_conv.dat > convergence_KP.dat
BEST_KPOINTS=\$(./convergence_summary kpoints.dat convergence_KP.dat )
echo  \$BEST_KPOINTS > BEST_KPOINTS.dat
rm __*.dat

hostname
date
echo "#######################################################"
echo "###  VASP calculation for KPOINTS convergence ends! ###"
echo "#######################################################"
EOF
echo "run_vasp_KP-conv.sh file has been successfully created"

