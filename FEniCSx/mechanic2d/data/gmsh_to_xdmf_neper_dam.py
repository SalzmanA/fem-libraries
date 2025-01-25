from mpi4py import MPI
from dolfinx.io import XDMFFile, gmshio
msh, c, e = gmshio.read_from_msh("../../../common/data/neper_dam.msh_2.2", MPI.COMM_WORLD, gdim=2)
msh.name = "neper_dam"
c.name = f"{msh.name}_cells"
e.name = f"{msh.name}_facets"
with XDMFFile(msh.comm, "neper_dam.xdmf", "w") as file:
    #msh.topology.create_connectivity(2, 3)
    file.write_mesh(msh)
    file.write_meshtags(
        c, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
    )
    file.write_meshtags(
        e, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
    )
