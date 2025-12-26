from ase.build import bulk
from ase import Atoms
from ase.io import write
import numpy as np


lattice_constant = 3.615  
cutoff_radius = 10.5  
simulation_box_size = 200.0  


fcc_bulk = bulk("Cu", "fcc", a=lattice_constant, cubic=True) * (50, 50, 50)


bulk_positions = fcc_bulk.get_positions()
bulk_center = bulk_positions.mean(axis=0)


positions = []
symbols = []
for atom in fcc_bulk:
    x, y, z = atom.position
    
    distance = np.sqrt((x - bulk_center[0])**2 + (y - bulk_center[1])**2 + (z - bulk_center[2])**2)
    
    if distance <= cutoff_radius:
        positions.append([x, y, z])
        symbols.append("Cu")


cluster = Atoms(symbols=symbols, positions=positions)


cluster.center(vacuum=simulation_box_size / 2 - cutoff_radius)


output_filename = "CNP_2nm.xyz"
write(output_filename, cluster)