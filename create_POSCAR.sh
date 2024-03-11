#!/bin/bash

# コマンドライン引数から.vaspファイルのパスを取得
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <VASP_FILE_PATH>"
    exit 1
fi
FILE_PATH="$1"

# .vaspファイルの存在を確認
if [ ! -f "$FILE_PATH" ]; then
    echo "Error: File does not exist - $FILE_PATH"
    exit 2
fi

# .vaspファイルをPOSCARにコピー
cp "$FILE_PATH" POSCAR

echo "POSCAR file has been successfully created from $FILE_PATH"

