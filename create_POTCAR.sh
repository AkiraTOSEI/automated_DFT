#!/bin/bash

# コマンドライン引数の処理
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <VASP_FILE_PATH> [precision]"
    exit 1
fi
FILE_PATH="$1"
PRECISION="${2:-high}" # precision引数がなければ、デフォルトで'high'を使用

# ファイルの存在を確認
if [ ! -f "$FILE_PATH" ]; then
    echo "Error: File does not exist - $FILE_PATH"
    exit 2
fi

# 6行目から原子のリストを取得する前にフラグを初期化
first_atom=true

# 6行目から原子のリストを取得
ATOMS=$(awk 'NR==6 {print $0}' "$FILE_PATH")

# ATOMSが空かどうかを確認
if [ -z "$ATOMS" ]; then
    echo "Error: No atoms found in the file - $FILE_PATH"
    exit 3
fi

# used_POTCAR.txtを初期化
echo "" > used_POTCAR.txt

# 原子ごとにループ
for ATOM in $ATOMS; do
    # precisionに応じた優先度リストを設定
    case $PRECISION in
        high)
            SUFFIXES=("${ATOM}_d" "${ATOM}_sv" "${ATOM}_pv" "$ATOM")
            ;;
        mid)
            SUFFIXES=("${ATOM}_d" "${ATOM}_pv" "$ATOM" "${ATOM}_sv")
            ;;
        low)
            SUFFIXES=("$ATOM" "${ATOM}_d" "${ATOM}_pv" "${ATOM}_sv")
            ;;
        *)
            echo "Error: Invalid precision - $PRECISION"
            exit 4
            ;;
    esac

    # 優先度に基づいてファイルを検索
    for SUFFIX in "${SUFFIXES[@]}"; do
        POTCAR_PATH="/home/share/VASP/potpaw_PBE.54/${SUFFIX}/POTCAR"
        if [ -f "$POTCAR_PATH" ]; then
            echo "$POTCAR_PATH" >> used_POTCAR.txt
            if [ "$first_atom" = true ]; then
                cat "$POTCAR_PATH" > POTCAR
                first_atom=false
            else
                cat "$POTCAR_PATH" >> POTCAR
            fi
            break
        fi
    done
done

echo "POTCAR file has been successfully created/updated."
echo "Used POTCAR paths have been recorded in used_POTCAR.txt."
