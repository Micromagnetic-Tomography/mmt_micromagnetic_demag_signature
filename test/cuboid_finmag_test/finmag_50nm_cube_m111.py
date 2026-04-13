import finmag
import dolfin as df
import numpy as np

x0, y0, z0 = -25, -25, -60
x1, y1, z1 = 25, 25, -10
box = finmag.util.mesh_templates.Box(x0, y0, z0, x1, y1, z1)

mesh = box.create_mesh(maxh=2.0, save_result=False)

# Continuous Galerkin so we have 1 value per node
V = df.FunctionSpace(mesh, 'CG', 1)
v = df.TestFunction(V)
tet_vol = df.assemble(v * df.dx)

# unit_length = 1e-9
# dim = 3
# vols = df.assemble(v * df.dx).array() * unit_length ** dim
# print(vols)
# print(np.sum(vols))

sim = finmag.Simulation(mesh, Ms=1e6, unit_length=1e-9)
sim.m = (1.0, 1.0, -1.0)

# sim.m is obtained in xxx format: x1 x2 ... y1 y2 ... z1 z2 ...
# coordinates and vols in nm units
data = np.column_stack((sim.mesh.coordinates(), sim.m.reshape(3, -1).T, tet_vol.array()))
np.savetxt('finmag_cube_m111.txt', data)
# print(data[:, -1].sum())
