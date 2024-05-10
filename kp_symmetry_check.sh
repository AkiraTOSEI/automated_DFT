#!/bin/bash

# Read data from BEST_KPOINTS.dat
read kx ky kz < ./BEST_KPOINTS.dat

# Check if the data contains three integers
if ! [[ "$kx" =~ ^[0-9]+$ && "$ky" =~ ^[0-9]+$ && "$kz" =~ ^[0-9]+$ ]]; then
  echo "Error: The file does not contain three integers."
  exit 1
fi

# Determine if all three numbers are equal
if [ "$kx" -eq "$ky" ] && [ "$ky" -eq "$kz" ]; then
  # If all numbers are equal, create an empty file
  > INCAR_tail
else
  # If numbers are not equal, write specified content to the file
  echo "ISYM = -1" > INCAR_tail
  echo "SYMPREC = 0.00000001" >> INCAR_tail
fi
