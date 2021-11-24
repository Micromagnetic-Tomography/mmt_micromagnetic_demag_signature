# Using modified/forked meshio:
# https://github.com/davidcortesortuno/meshio
import meshio
# import numpy as np
# import re
import sys
from pathlib import Path


FNAME = Path(sys.argv[1])

gmsh_mesh = meshio.read(FNAME)

mesh = meshio.Mesh(
            gmsh_mesh.points, gmsh_mesh.cells,
            cell_data=dict(SD=gmsh_mesh.cell_data['gmsh:geometrical'])
            )

# Notice that MERRILL 1.3.2 expects different headers so we must manually
# change: DATAPACKING -> F , ZONETYPE -> ET , BLOCK -> FEBLOCK
#         FETETRAHEDRON -> TETRAHEDRON
meshio.write(FNAME.parents[0] / Path(FNAME.stem + '_gmsh_NOMAG.tec'), mesh,
             data_formats=dict(SD='{:.0f}'))
# meshio.write(FNAME.parents[0] / Path(FNAME.stem + '_gmsh.xml'), mesh)
