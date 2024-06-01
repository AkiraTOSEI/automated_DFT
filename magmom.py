def process_poscar(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    # Read the 6th line (index 5)
    atoms_line = lines[5].strip().split()
    num_line = lines[6].strip().split()
    
    print("MAGMOM =", end="")
    
    # Process each atom
    for atom, num  in zip(atoms_line, num_line):
        if atom in ['Fe', 'Co', 'Ni', 'Gd']:
            mag_val = 1.0
        else:
            mag_val = 0.0
        for _ in range(int(num)):
            print(f" {mag_val}", end="")
    
    print("")

process_poscar('POSCAR_gt')
