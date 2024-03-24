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


# 収束判定を初期化
result=FALSE

EOF

i=1
while IFS=, read -r kx ky kz; do

    cat >> run_vasp_KP-conv.sh <<EOF
if [ "\$result" == "True" ]; then
    echo "KPOINT Converged."
else
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
    # 何らかのエラーでエネルギーが取得できなかった場合は0を入れる。
    if [[ -z "\${energy// }" ]]; then
        energy="0.0" 
    fi 
fi
EOF

    if [ $i -eq 1 ]; then
        cat >> run_vasp_KP-conv.sh <<EOF
echo "$kx $ky $kz \$energy" > kptest.dat
echo "\$energy" > kp_e.dat
EOF
    else
        cat >> run_vasp_KP-conv.sh <<EOF
if [ "\$result" == "True" ]; then
    echo "KPOINT Converged."
else
    echo "$kx $ky $kz \$energy" >> kptest.dat
    echo "\$energy" >> kp_e.dat
    result=\$(./convergence_checker_KP kp_e.dat)
fi
EOF
    fi
    i=$((i+1))

done < kpoints_candidate.csv

cat >> run_vasp_KP-conv.sh <<EOF
hostname
date
echo "#######################################################"
echo "###  VASP calculation for KPOINTS convergence ends! ###"
echo "#######################################################"
EOF
echo "run_vasp_KP-conv.sh file has been successfully created"

