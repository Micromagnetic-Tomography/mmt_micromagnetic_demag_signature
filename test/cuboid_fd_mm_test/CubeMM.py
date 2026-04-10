import numpy as np

from mumaxplus import Ferromagnet, Grid, World
from mumaxplus.util import plot_field
import mumaxplus.util.shape as shapes
from mumaxplus.util.show import show_field_3D

nm = 1e-9
mu0 = 4 * np.pi * 1e-7
# Material parameters
# A = 1e-12
# Ms = 1e6

# system dimensions
lxy = 50.0 * nm
dxy = 2 * nm
dz = 2 * nm
lz = 50.0 * nm
nxy = int(lxy / dxy)
nz = int (lz / dz)

# geometry parameters
cube = shapes.Cube(lxy).translate(lxy/2, lxy/2, lz/2)

# Create the world
grid_size = (nxy, nxy, nz)
cell_size = (dxy, dxy, dz)
world = World(cell_size)

# Ferromagnet
magnet = Ferromagnet(world, Grid(size=grid_size), geometry=cube, name="cube_mm")
magnet.magnetization = (1.0, 1.0, -1.0)

# plot_field(magnet.magnetization, arrow_size=8, layer=0)

magnet_data = np.array(magnet.magnetization.eval(), dtype=float)
mesh_data = np.array(magnet.magnetization.meshgrid, dtype=float)  # (3, Nz, Ny, Nx)
magnet_full = np.concatenate([mesh_data.reshape(3, -1).T, 
                              magnet_data.reshape(3, -1).T], axis=1)
# Scale to nm
magnet_full[:, :3] *= 1e9

# manually translate z in -60 nm
magnet_full[:, 0] -= (25.0 - dxy * 1e9 * 0.5)
magnet_full[:, 1] -= (25.0 - dxy * 1e9 * 0.5)
magnet_full[:, 2] -= (60.0 - dz * 1e9 * 0.5)

name_m = "cube_mumaxPlus_L_50nm_centerAt_-35nm_dxyz_2nm.npy"
np.save(name_m, magnet_full)
