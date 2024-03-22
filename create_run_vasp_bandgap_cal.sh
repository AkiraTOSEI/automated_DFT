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
cp POTCAR KPOINTS calculate_bandgap.cpp ./bandgap_cal
cd bandgap_cal

echo "################################################"
echo "################################################"
echo "###  VASP calculation for  Bandgap starts!  ###"
echo "################################################"
echo "################################################"
echo ""

# ENCUTの収束後値をもとに、bandgap計算用のINCARを作成する.
ENCUT=\$(tail -n -1 ../encut_cutoff.dat | awk '{print \$1}')
echo "bandap calculation INCAR" > INCAR
echo "ISTART = 1" >> INCAR
echo "ICHARG = 11" >> INCAR
echo "ISPIN = 1" >> INCAR
echo "ENCUT = \$ENCUT" >> INCAR
echo "EDIFF = 1.0e-4" >> INCAR
echo "LWAVE = .F." >> INCAR
echo "ISMEAR = -5" >> INCAR
echo "#SIGMA = 0.01" >> INCAR
echo "PREC = accurate" >> INCAR
echo "LREAL = .F." >> INCAR
echo "NELM = 200" >> INCAR
echo "NPAR = 4" >> INCAR
echo "LORBIT = 11" >> INCAR
echo "NEDOS = 2000" >> INCAR

mpiexec -iface ib0 -launcher rsh -machinefile \$PBS_NODEFILE -ppn 16 /home/share/VASP/vasp.5.4.4

echo "#######################################"
echo "###            VASP ends!           ###"
echo "#######################################"

cp OUTCAR ../final_OUTCAR
cp DOSCAR ../final_DOSCAR
cp POSCAR ../final_POSCAR
cp KPOINTS ../final_KPOINTS
cp INCAR ../final_INCAR

g++ -o calculate_bandgap calculate_bandgap.cpp -std=c++11
E_fermi=\$(grep "E-fermi" ../final_OUTCAR | awk '{print \$3}')
OUTPUT_FILE="../../bandgap_result.csv"

echo "E_fermi: \$E_fermi"

# outputが存在するかどうかをチェック
if [ -e \$OUTPUT_FILE ]; then
    # ファイルが存在する場合、通常のlsの結果をファイルに出力
    ./calculate_bandgap \$E_fermi >> \$OUTPUT_FILE
else
    # ファイルが存在しない場合、隠しファイルを含むls -aの結果をファイルに出力
    echo 'file_name, bandgap, nonmetal' > \$OUTPUT_FILE
    ./calculate_bandgap \$E_fermi >> \$OUTPUT_FILE
fi

echo "#######################################"
echo "###     bandgap cpp ends!           ###"
echo "#######################################"

cd ../

# 特定のファイル以外を削除
rm -rf bandgap_cal final_sc
shopt -s extglob
rm -f !(final_INCAR|final_KPOINTS|final_OUTCAR|final_DOSCAR|final_POSCAR|used_POTCAR.txt|kptest.dat|encut_test.dat)

date



EOF
echo "run_vasp_bandgap.sh file has been successfully created"





# bandgap計算ファイルの内容を定義
cat > calculate_bandgap.cpp << EOF
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "使用法: " << argv[0] << " eFermi" << std::endl;
        return 1;
    }   
    
    double eFermi = std::stod(argv[1]); // コマンドラインからeFermiの値を受け取る
    std::string doscarPath = "../final_DOSCAR"; // DOSCARファイルのパスを適宜設定してください
    std::ifstream doscarFile(doscarPath);
    std::string line;
    double a = 0, b = 0, c = 0;
    bool nonmetal = false;

    if (!doscarFile.is_open()) {
        std::cerr << "ファイルを開けませんでした: " << doscarPath << std::endl;
        return 1;
    }   

    double closestBelowValue = -std::numeric_limits<double>::max();
    double smallestAboveZeroSecondCol = std::numeric_limits<double>::max();
    double closestAboveValue = std::numeric_limits<double>::max();

    while (getline(doscarFile, line)) {
        std::istringstream iss(line);
        double firstCol, secondCol;
        if (!(iss >> firstCol >> secondCol)) { continue; }

        if (firstCol < eFermi && firstCol > closestBelowValue) {
            closestBelowValue = firstCol;
            a = firstCol;
        }   
        if (firstCol > eFermi && secondCol > 0 && firstCol < smallestAboveZeroSecondCol) {
            smallestAboveZeroSecondCol = firstCol;
            b = firstCol;
        }   
        if (firstCol > eFermi && firstCol < closestAboveValue) {
            closestAboveValue = firstCol;
            c = secondCol;
        }   
    }   

    doscarFile.close();

    double bandgap = b - a;
    nonmetal = c == 0;
                                                                                                                                                                                                                                                                                                                                
    std::cout << "$FILE_PATH , " << bandgap << ", " << nonmetal << std::endl;

    return 0;
}
EOF


echo "calculate_bandgap.cpp file has been successfully created"
