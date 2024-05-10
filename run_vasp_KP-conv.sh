#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N POSCAR
cd $PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=$(cat $PBS_NODEFILE | wc -l)




echo "########################################################"
echo "###  VASP calculation for KPOINTS convergence stars! ###"
echo "########################################################"
echo "########################################################"

cp POSCAR_distorted POSCAR

# 収束判定を初期化
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 4 4 4
echo "Processing k-point 4 4 4"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "4 4 4" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "4 4 4"  > kpoints.dat
echo "$energy" > energy_KP.dat
echo $FORCE > force_KP.dat 
echo $PRESSURE > pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 5 5 5
echo "Processing k-point 5 5 5"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "5 5 5" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "5 5 5" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 6 6 6
echo "Processing k-point 6 6 6"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "6 6 6" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "6 6 6" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 7 7 7
echo "Processing k-point 7 7 7"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "7 7 7" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "7 7 7" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 8 8 8
echo "Processing k-point 8 8 8"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "8 8 8" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "8 8 8" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 9 9 9
echo "Processing k-point 9 9 9"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "9 9 9" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "9 9 9" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 10 10 10
echo "Processing k-point 10 10 10"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "10 10 10" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "10 10 10" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 11 11 11
echo "Processing k-point 11 11 11"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "11 11 11" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "11 11 11" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 12 12 12
echo "Processing k-point 12 12 12"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "12 12 12" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "12 12 12" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 13 13 13
echo "Processing k-point 13 13 13"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "13 13 13" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "13 13 13" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 14 14 14
echo "Processing k-point 14 14 14"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "14 14 14" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "14 14 14" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 15 15 15
echo "Processing k-point 15 15 15"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "15 15 15" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "15 15 15" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 16 16 16
echo "Processing k-point 16 16 16"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "16 16 16" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "16 16 16" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat
# Append k-points processing commands to run_vasp_KP-conv.sh
### create KPOINTS file
# K-point: 17 17 17
echo "Processing k-point 17 17 17"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "17 17 17" >> KPOINTS
echo "0 0 0" >> KPOINTS

### calculcation
mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4
energy=$(grep 'free  ' OUTCAR | tail -1 | awk '{print $5}')
FORCE=$(awk '/TOTAL-FORCE/ {getline; getline; print $4, $5, $6}' OUTCAR) 
PRESSURE=$(grep "external pressure" OUTCAR | awk '{print $4}')
# 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
if [[ -z "${energy// }" ]]; then
    energy="0.0" 
else
    cp KPOINTS final_KPOINTS
fi

if [[ -z "${FORCE// }" ]]; then
    FORCE="1.0 1.0 1.0" 
fi


if [[ -z "${PRESSURE// }" ]]; then
    PRESSURE="100.0" 
fi

echo "17 17 17" >> kpoints.dat
echo "$energy" >> energy_KP.dat
echo $FORCE >> force_KP.dat 
echo $PRESSURE >> pressure_KP.dat


#収束判定ができるcppファイルをコンパイルする 
g++ -std=c++11 -o check_conv convergence_checker.cpp
g++ -std=c++11 -o calculate_l2_norm l2_norm.cpp
g++ -std=c++11 -o calculate_vol calculate_vol.cpp
g++ -std=c++11 -o convergence_summary convergence_summary.cpp

# 原子数を取得する
NUM_ATOM=$(sed -n '7p' POSCAR | awk '{sum=0; for (i=1; i<=NF; i++) sum+=$i; print sum}')

# Pressureの最終値との差分が5 GPa, 絶対値が5 GPa以下になるようにする. 
# VASPでは圧力単位がKilo Barなので、GPaにへんかんするために0.1をかける
./check_conv pressure_KP.dat 0.1 0.5 last > pressure_diff.dat
./check_conv pressure_KP.dat 0.1 5. zero > pressure_abs.dat
# エネルギー/原子の最終値との差分が1meV以下になるようにする
./check_conv energy_KP.dat "1000/$NUM_ATOM" 1 last > energy_KP_diff.dat
# 力はl2_normをとり、その値が
./calculate_l2_norm force_KP.dat > force_norm.dat
./check_conv force_norm.dat 1. 0.01 last > force_diff.dat

# KPOINTの体積を計算する.
./calculate_vol kpoints.dat > k_vol.dat

# 結果ファイルを作成する
paste force_diff.dat pressure_abs.dat pressure_diff.dat energy_KP_diff.dat > __tmp_conv.dat
echo "Fx Fy Fz F_norm(eV/Å) F_diff(eV/Å) pressure(GPa) P_diff(GPa) Energy(meV/atom) E_diff(meV/atom)" > __kp_history.dat
awk '{print $1, $2, $7, $8, $10, $11}' __tmp_conv.dat > tmp.dat; paste force_KP.dat tmp.dat| expand | tr -s "[:space:]" >> __kp_history.dat
echo  "kx ky kz k_vol" > kp_tmp1.dat
paste kpoints.dat k_vol.dat > kp_tmp2.dat 
cat kp_tmp1.dat kp_tmp2.dat > kp_tmp.dat
paste kp_tmp.dat __kp_history.dat | expand | tr -s "[:space:]" >  kp_history.dat
# 収束判定のファイルを作成する
awk '{print $3, $6, $9, $12}' __tmp_conv.dat > convergence_KP.dat
BEST_KPOINTS=$(./convergence_summary kpoints.dat convergence_KP.dat )
echo  $BEST_KPOINTS > BEST_KPOINTS.dat
rm __*.dat

# BEST_KPOINTSファイルの存在、内容、及び数値のチェック
if [ ! -f BEST_KPOINTS.dat ] || [ ! -s BEST_KPOINTS.dat ] || ! grep -q '[0-9]' BEST_KPOINTS.dat; then
  echo "エラー: 'BEST_KPOINTS.dat' ファイルが存在しないか、中身が空、または数値が含まれていません。"
  exit 1
else
  echo "'BEST_KPOINTS.dat' ファイルは正常に存在し、中身に数値が含まれています。"
  # ここにBEST_KPOINTSファイル用の処理を書く
fi


### BEST_KPOINTSの対称性によって、INCARに追記すべき項目を作成する
sh kp_symmetry_check.sh


hostname
date
echo "#######################################################"
echo "###  VASP calculation for KPOINTS convergence ends! ###"
echo "#######################################################"
