#!/bin/bash








##########################################
# OUTCAR, DOSCAR からbandgapを計算する ###
#########################################


# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <VASP_FILE_PATH>"
    exit 1
fi
FILE_PATH="$1"



# INCARファイルの内容を定義
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
    std::string doscarPath = "final_DOSCAR"; // DOSCARファイルのパスを適宜設定してください
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


g++ -o calculate_bandgap calculate_bandgap.cpp -std=c++11
E_fermi=$(grep "E-fermi" ./bandgap_cal/OUTCAR | awk '{print $3}')
OUTPUT_FILE='../bandgap_result.csv'

# outputが存在するかどうかをチェック
if [ -e $OUTPUT_FILE ]; then
    # ファイルが存在する場合、通常のlsの結果をファイルに出力
    ./calculate_bandgap $E_fermi >> $OUTPUT_FILE
else
    # ファイルが存在しない場合、隠しファイルを含むls -aの結果をファイルに出力
    echo 'file_name, bandgap, nonmetal' > $OUTPUT_FILE
    ./calculate_bandgap $E_fermi >> $OUTPUT_FILE
fi


