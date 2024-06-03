#!/bin/bash




# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <VASP_FILE_PATH>"
    exit 1
fi
FILE_PATH="$1"


cat > run_vasp_bandgap.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd \$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=\$(cat \$PBS_NODEFILE | wc -l)



#################################
######     bandgap       ##########
#################################


# 最終計算用のディレクトリを作り、必要なファイルを移動させる。
mkdir -p bandgap_cal
cp POTCAR calculate_bandgap.py ./bandgap_cal
cd bandgap_cal

echo "################################################"
echo "################################################"
echo "###  VASP calculation for  Bandgap starts!  ###"
echo "################################################"
echo "################################################"
echo ""

# ENCUTの収束後値をもとに、sc計算用のINCARを作成する.
BEST_KPOINTS=\$(cat ../BEST_KPOINTS.dat)
MESH_METHOD=\$(cat ../mesh_method.dat)
echo "BEST KPOINTS"
echo "k-points" > KPOINTS
echo "0" >> KPOINTS
echo "\$MESH_METHOD" >> KPOINTS
echo "\$BEST_KPOINTS" >> KPOINTS
echo "0 0 0" >> KPOINTS


# ENCUTの収束後値をもとに、bandgap計算用のINCARを作成する.
BEST_ENCUT=\$(cat ../BEST_ENCUT.dat)
EDIFF=\$(cat ../EDIFF.dat)
echo "bandap calculation INCAR" > INCAR
echo "ISTART = 1 #default:1, 0:scratch, 1:WAVECARを読む、なければ0になる " >> INCAR
echo "ICHARG = 11 #default:2(ISTART=0)0,else:0, 0:初期波動関数から電荷を計算,2:原子電荷の重ね合わせを使う,11:CHGCAR利用" >> INCAR
echo "ISPIN = 2 # default:1, 1:スピン分極なし、2:スピン分極あり" >> INCAR
echo "\$BEST_ENCUT" >> INCAR
echo "LWAVE = .F. # WAVECARを出力するか" >> INCAR
echo "ISMEAR = -5 # default:1,部分占有軌道の設定、絶縁体・半導体では-5を推奨 " >> INCAR
echo "#SIGMA = 0.01 # default:0.2, smearingの幅、系の事前知識がないときはISMEAR=0でSIGMA=0.03~0.05 " >> INCAR
echo "PREC = accurate" >> INCAR
echo "LREAL = .F. # 精度の指定。default NORMAL" >> INCAR
echo "NELM = 200 # default:60, SC計算の最大回数を指定" >> INCAR
#echo "NPAR = 4 # defalut:コア数、並列計算に使うコア数" >> INCAR
echo "LORBIT = 11 # default: None, 11:DOSCAR, PROCARを出力" >> INCAR
echo "NEDOS = 2000 #default:301, DOS計算のグリッド数" >> INCAR
echo \$EDIFF >> INCAR
cat ../INCAR_tail >> INCAR
cat ../incar_magmom.dat >> INCAR

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4

echo "#######################################"
echo "###            VASP ends!           ###"
echo "#######################################"

cp OUTCAR ../final_OUTCAR
cp DOSCAR ../final_DOSCAR
cp POSCAR ../final_POSCAR
cp KPOINTS ../final_KPOINTS
cp INCAR ../final_INCAR

### バンドギャップ計算用のDOSCARを作成する
# 最初の8行をスキップし、4列以上のデータを含む行を無視
tail -n +9 DOSCAR | awk 'NF <= 5 { print }' > tmp_DOSCAR

E_fermi=\$(grep "E-fermi" ../final_OUTCAR | awk '{print \$3}')
OUTPUT_FILE="../../bandgap_result.csv"

echo "E_fermi: \$E_fermi"

### E-fermiの値が計算できたかどうかでダミーの値をいれるかどうかを判定する
if [ -z "\$E_fermi"]; then
    # outputが存在するかどうかをチェック
    if [ -e \$OUTPUT_FILE ]; then
        # ファイルが存在する場合、通常のlsの結果をファイルに出力
        echo "$FILE_PATH , -1, False" >> \$OUTPUT_FILE
    else
        # ファイルが存在しない場合、隠しファイルを含むls -aの結果をファイルに出力
        echo 'file_name, bandgap, nonmetal' > \$OUTPUT_FILE
        echo "$FILE_PATH , -1, False" >> \$OUTPUT_FILE
    fi
else
    # outputが存在するかどうかをチェック
    if [ -e \$OUTPUT_FILE ]; then
        # ファイルが存在する場合、通常のlsの結果をファイルに出力
        python calculate_bandgap.py -f \$E_fermi >> \$OUTPUT_FILE
    else
        # ファイルが存在しない場合、隠しファイルを含むls -aの結果をファイルに出力
        echo 'file_name, bandgap, nonmetal' > \$OUTPUT_FILE
        python calculate_bandgap.py -f \$E_fermi >> \$OUTPUT_FILE
    fi
fi


echo "#######################################"
echo "###     bandgap cpp ends!           ###"
echo "#######################################"

cd ../

# 特定のファイル以外を削除
rm -rf bandgap_cal final_sc relax_cal
shopt -s extglob
rm -f !(POSCAR_gt|POSCAR_distorted|final_INCAR|final_KPOINTS|final_OUTCAR|final_DOSCAR|final_POSCAR|used_POTCAR.txt|kp_history.dat|encut_history.dat)

date



EOF
echo "run_vasp_bandgap.sh file has been successfully created"





# bandgap計算ファイルの内容を定義
cat > calculate_bandgap.py << EOF
import pandas as pd
import argparse

def analyze_dos(e_fermi=None):
    if e_fermi is None:
        raise ValueError("Fermi energy must be specified.")
    
    # Load data
    dos = pd.read_csv('tmp_DOSCAR', delim_whitespace=True, header=None)
    dos = dos.rename(columns={0: 'Energy', 1: 'DOS(UP)', 2: 'DOS(DOWN)', 3: 'Integ_DOS(UP)', 4: 'Integ_DOS(DOWN)'})
    
    # Calculate if energy values are above the Fermi level
    dos['AboveEf'] = dos['Energy'] - float(e_fermi) > 0
    
    # Get the energy, DOS(UP), DOS(DOWN) of the first point above the Fermi level
    energy, dup, ddown = dos[dos['AboveEf']].head(1)[['Energy', 'DOS(UP)', 'DOS(DOWN)']].values.tolist()[0]
    
    if dup > 0 or ddown > 0:
        non_metal = False
        band_gap = 0
    else:
        non_metal = True
        
        # Get Energy of the top of the valence band
        top_VB_idx = dos[dos['AboveEf']].index[0] - 1
        top_VB_energy = dos['Energy'].loc[top_VB_idx]
        
        # Get Energy of the bottom of the conduction band
        dos_above_Ef = dos[dos['AboveEf']][['Energy', 'DOS(UP)', 'DOS(DOWN)']]
        dos_above_Ef['state_exist'] = (dos_above_Ef[['DOS(UP)', 'DOS(DOWN)']] > 0).any(axis=1)
        bottom_CB_energy = dos_above_Ef[dos_above_Ef['state_exist']].head(1)['Energy'].values[0]
        
        band_gap = bottom_CB_energy - top_VB_energy
    
    print('$FILE_PATH' ,f', {band_gap:.3f}, {non_metal}')

def main():
    parser = argparse.ArgumentParser(description="Analyze DOS data for band gaps.")
    parser.add_argument("-f", "--fermi", type=float, required=True, help="Fermi energy level")
    
    args = parser.parse_args()
    analyze_dos(args.fermi)

if __name__ == "__main__":
    main()
EOF

echo "calculate_bandgap.py file has been successfully created"
