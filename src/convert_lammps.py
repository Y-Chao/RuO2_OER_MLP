import sys

from ase.data import atomic_numbers
from ase.io import lammpsdata, read

# def out_fixed_index(index, c_index):

print("Usage: python convert_lammps.py <input_xyz> <output_lammps_data>")

struc = read(sys.argv[1])
lmp_file = sys.argv[2]

order = list(set(struc.symbols))
order.sort(key=lambda x: atomic_numbers[x])

lammpsdata.write_lammps_data(lmp_file, struc, specorder=order)

atoms_index = [i.index for i in struc]
try:
    c_index = struc.constraints[0].index
    uc_index = [str(i + 1) for i in atoms_index if i not in c_index]
    cc_index = [str(i + 1) for i in c_index]

    print("Constraint index: ")
    cc_sym = " ".join(cc_index)
    print(cc_sym)

    print("Unonstraint index: ")
    uc_sym = " ".join(uc_index)
    print(uc_sym)
except:
    print("No constraint")
