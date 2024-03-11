#!/bin/bash

# INCARファイルのパスを指定
INCAR_FILE="INCAR"

# INCARファイルの内容を定義
cat > $INCAR_FILE << EOF
System = default
ISTART = 0 ; ICHARG = 2
ENCUT = 240
ISMEAR = -5; SIGMA = 0.1
EOF

echo "INCAR file has been successfully created."

