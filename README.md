# 概要
POSCARを入力とし、以下の手順を自動化するスクリプト群。qsubするスクリプト内で以下の計算を全て記述しないといけないため、「スクリプトを書くスクリプト」でスクリプトを書いたあとに、連結する構成にしている。分割したのは、デバッグのため。
1. K点の収束計算
2. INCARにおけるENCUTの収束
3. 緩和計算
4. バンドギャップの計算

# 使い方
`all_convergence.sh`に引数として、POSCARのパスを入れる。e.g) `sh all_convergence.sh poscar_file/normal_Si`

# ファイルの内容
- `all_convergence.sh`: 計算を回すためのファイル
- `create_INCAR.sh` : K点収束計算のためのINCARファイルを生成するスクリプト
- `create_POTCAR.sh`: `all_convergence.sh`の引数として入れたPOSCARファイルから、POTCARファイルを構成するためのスクリプト
- `create_POSCAR.sh`: `all_convergence.sh`の引数として入れたPOSCARファイルを、ファイル名をPOSCARに変更するためのスクリプト
- `create_conv_cpp.sh` : 収束判定を行うためのC++ファイルを生成するためのスクリプト
- `generate_kpoints.py` : `all_convergence.sh`の引数として入れたPOSCARファイルから、結晶ベクトルの比をもとにK点収束用のK点リストを生成するpythonファイル
- `create_run_vasp_KP-conv.sh` : K点収束計算を行うスクリプトを生成するスクリプト
- `create_run_vasp_ENCUT-conv.sh` : INCARにおけるENCUTの収束計算を行うスクリプトを生成するスクリプト
- `create_run_vasp_relax.sh` : 緩和計算を行うスクリプトを生成するスクリプト
- `create_run_vasp_bandgap_cal.sh` : バンドギャップ計算のためのDFTを行うスクリプトを生成するスクリプト
- `calculate_bandgap.sh` : バンドギャップ計算を行ったDOSCARから、バンドギャップを計算するC＋＋ファイルを生成し、バンドギャップ計算を行うスクリプト

# 結果ファイル
- `bandgap_result.csv` : bandgapの計算結果が入ったファイル。なお、nonmetal列が1ならば、bandgapはゼロである。（入力された値はDOSCARにおけるステップ幅なので、無意味）。例として入っているfcc_Siは導体なので、nonmetal列が1になっている。
- `workspace__*/kptest.dat` : K点収束計算の結果
- `workspace__*/encut_test.dat` :  INCARにおけるENCUTの収束計算の結果
- `workspace__*/final_DOSCAR` : 緩和計算後の最終的なバンドギャップ計算を行った後のDOSCARファイル。bandgap_result.csvは、このファイルをもとに計算する
- `workspace__*/final_OUTCAR` : 緩和計算後の最終的なバンドギャップ計算を行った後のOUTCARファイル。
- `workspace__*/used_POTCAR.txt` : どのPOTCARファイルを使ったかを記録するためのファイル。
  
