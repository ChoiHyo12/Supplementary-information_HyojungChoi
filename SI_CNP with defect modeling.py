from ase.build import bulk
from ase import Atoms
from ase.io import write
import numpy as np
import random


lattice_constant = 3.615  
cutoff_radius = 20.0  
simulation_box_size = 200.0  


fcc_bulk = bulk("Cu", "fcc", a=lattice_constant, cubic=True) * (30, 30, 30)


bulk_center = fcc_bulk.get_positions().mean(axis=0)


positions = []
symbols = []
for atom in fcc_bulk:
    x, y, z = atom.position
    
    distance = ((x - bulk_center[0])**2 + (y - bulk_center[1])**2 + (z - bulk_center[2])**2)**0.5
    
    if distance <= cutoff_radius:
        positions.append([x, y, z])
        symbols.append("Cu")


cluster = Atoms(symbols=symbols, positions=positions)


cluster.center(vacuum=simulation_box_size / 2 - cutoff_radius)


cluster_center = cluster.get_positions().mean(axis=0)
positions = cluster.get_positions()


distances = np.linalg.norm(positions - cluster_center, axis=1)


surface_threshold = 0.85 * cutoff_radius  
surface_mask = distances >= surface_threshold
internal_mask = distances < surface_threshold


internal_indices = np.where(internal_mask)[0]
surface_indices = np.where(surface_mask)[0]


num_to_remove = int(0.02 * len(internal_indices)) 
indices_to_remove = random.sample(list(internal_indices), num_to_remove)


remaining_indices = np.setdiff1d(range(len(cluster)), indices_to_remove)


defected_cluster = cluster[remaining_indices]


output_filename = "defect4_2%.xyz"
write(output_filename, defected_cluster)

