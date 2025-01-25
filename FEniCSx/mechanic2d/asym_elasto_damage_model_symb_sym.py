#import pdb
from scipy import sparse
from mpi4py import MPI
start_all=MPI.Wtime()
start=MPI.Wtime()
import numpy as np
from pathlib import Path
from ctypes import CDLL
import sympy as sp
import ufl
#from ufl.classes import Sqrt as sqrt
from ufl.classes import (Mesh,FunctionSpace,Coefficient,Constant,
                         TrialFunction,TestFunction,FacetNormal)
from ufl import ( sqrt,grad,inner,tr,Identity,dot,dx,ds,
                 conditional,lt,gt,eq,measure)
import basix
#from basix.ufl import element

from dolfinx.io import XDMFFile
from dolfinx.io import VTXWriter
from dolfinx import fem
from dolfinx import mesh
from dolfinx import la
from dolfinx import nls
from dolfinx import cpp
from dolfinx.fem import petsc
from dolfinx.nls import petsc
import adios4dolfinx
from petsc4py import PETSc

MAX_DAM=1.
MAX_REFINE=0
trac=True
dt=np.zeros(16)
do_print = (MPI.COMM_WORLD.Get_rank()==0)
nb_proc = MPI.COMM_WORLD.Get_size()
dim0=0
dim1=1
dim=2
dt[1]=MPI.Wtime()-start # 1. Initialize

# Read mesh =======================
start=MPI.Wtime()
#filename = Path("data/square.xdmf")
filename = Path("data/neper_dam.xdmf")
with XDMFFile(MPI.COMM_WORLD, filename, "r") as file:
    domain=file.read_mesh(ghost_mode=mesh.GhostMode.none,name='neper_dam')
    cell_tag=file.read_meshtags(domain,"neper_dam_cells")
    domain.topology.create_entities(dim1)
    edge_tag=file.read_meshtags(domain,"neper_dam_facets")
dt[2]=MPI.Wtime()-start # 2.1 Read the mesh

start=MPI.Wtime()
# refine ==========================
for i in range(MAX_REFINE): 
    ndo,parent_cell,parent_facet=mesh.refine_plaza(domain,edges=None,redistribute=False,option=cpp.refinement.RefinementOption.parent_cell_and_facet)
    cell_tag=mesh.transfer_meshtag(cell_tag,ndo,parent_cell)
    ndo.topology.create_entities(dim1)
    edge_tag=mesh.transfer_meshtag(edge_tag,ndo,parent_cell,parent_facet)
    domain=ndo
dt[3]+=MPI.Wtime()-start # 2.2 Refine the mesh

start = MPI.Wtime()
# finalize mesh topo ==============
topo=domain.topology
topo.create_entities(dim0) # for BC
topo.create_connectivity(dim0,dim)
topo.create_connectivity(dim1, dim0)
topo.create_connectivity(dim0, dim1)
edges=topo.connectivity(dim1,dim0)
dt[2]+=MPI.Wtime()-start # 2.1 Read the mesh

start=MPI.Wtime()
# element =========================
elem = basix.ufl.element("Lagrange", "triangle", 1, shape=(2, ))
element_scal = basix.ufl.element("Lagrange", "triangle", 1)
DGelement = basix.ufl.element("DG", "triangle", 0)
DGVelement = basix.ufl.element("DG", "triangle", 0, shape=(3,))
# spaces ==========================
space = fem.functionspace(domain, elem)
DGspace = fem.functionspace(domain, DGelement)
DGVspace = fem.functionspace(domain, DGVelement)
space_scal = fem.functionspace(domain, element_scal)
dt[5]+=MPI.Wtime()-start # 3.1 Define space

start=MPI.Wtime()
# damage ==========================
# damage field
d=fem.Function(space_scal,name="d")
d.x.array.fill(0.)
#dstart=MPI.Wtime()
# For ghost needs to put in place a correspondence of index between nodes and dofs as we
# gone use dofs/nodes for computation. To simplify a full correspondence is set for local
# and ghost. Pass by global looks not possible as python api do not show  global_indices.
# thus pass via ghosts which is a priori cheaper
# anyway it work only because Lagrange1 scalar dof and nodes are linked appropriately
im_space = space_scal.dofmap.index_map
im_mesh = topo.index_map(dim0)
nloc = space_scal.dofmap.index_map.size_local
nglob = space_scal.dofmap.index_map.num_ghosts
nlg = nloc+nglob
assert nloc == im_mesh.size_local,"Wrong node topo local"
assert nglob == im_mesh.num_ghosts,"Wrong node topo ghost"
corresp=np.arange(nlg,dtype=np.int32)
corresp[nloc:]=im_space.global_to_local(im_mesh.ghosts)
# glob not needed anymore 
# extra distributed vector to operate directly on mesh index
nd_valv=la.vector(im_mesh, bs=1)
nd_val=nd_valv.array
# number of edge per owned node
num_edges_per_nodes=np.zeros(nloc)
# collect edges chosen arbitrarily
#tag_edges_damaged = [4]
tag_edges_damaged = [148, 342, 333, 19,  380, 408, 328, 329, 325, 323, 96,  97,  531, 4,   471, 234, 235, 184, 236, 419, 350, 332, 364, 176, 77,  333, 341, 343, 144, 143]
edges_damaged=set()
for  tag in  tag_edges_damaged:
     s = edge_tag.find(tag)
     edges_damaged.update(s)
# collect nodes with damage
node_damaged=set()
for l in edges_damaged:
    link=edges.links(l)
    for i in link:
        node_damaged.add(corresp[i])
# set damage to MAX_DAM for node collected
d_val=d.x.array
for n in node_damaged:
    d_val[n]=MAX_DAM
# sum damage in owners
d.x.scatter_reverse(la.InsertMode.add)
# an owner may have damage above MAX_DAM if edges related to this owner are present in
# more than one process: reset do MAX_DAM
with np.nditer(d_val,op_flags=['readwrite']) as it:
    for di in it:
        if (di > MAX_DAM):
            di[...] = MAX_DAM
# propagate to ghost: a ghost may still be null if no edges are marked in current process but in other process ghost or
# owner counterparts are not null, thus this ghost has to be updated
d.x.scatter_forward()
# set edge counter
nodes=topo.connectivity(dim0,dim1)
nloc_edges = topo.index_map(dim1).size_local
rnlg=range(nlg)
for l in rnlg:
    link = nodes.links(l)
    nd_val[l] = 0.
    for ll in link:
      # edge owned ? yes, count contribution
      if (ll < nloc_edges):
          nd_val[l]+=1.
nd_valv.scatter_reverse(la.InsertMode.add)
rnloc=range(nloc)
for l in rnloc:
    if (nd_val[l]>0):
        num_edges_per_nodes[l] = 1./nd_val[l]
    else:
        num_edges_per_nodes[l] = 0.
#print("preamble")
#print(MPI.Wtime()-dstart)
#dstart=MPI.Wtime()
opRow=np.empty(nodes.array.size,dtype=np.int32)
opCol=np.empty(nodes.array.size,dtype=np.int32)
k=0
for l in rnlg:
    for ll in filter(lambda lf: lf<nloc_edges,nodes.links(l)):
        nn = edges.links(ll)
        vid = nn[0]
        if (vid == l): vid = nn[1]
        opCol[k]=corresp[vid]
        opRow[k]=l
        k=k+1
opRow=opRow[:k]
opCol=opCol[:k]
mat=sparse.csr_matrix((np.ones(k),(opRow,opCol)),shape=(nlg,nlg))
mat.sort_indices()
#print("matrix op",MPI.Wtime()-dstart)
# smooth damage around initial values
for iter_smooth in range(8*(MAX_REFINE+1)):
    nd_val.fill(0.)
    #dstart=MPI.Wtime()
    # first prod to enlarge
    nd_val[:]=mat.dot(d_val)
    toz=list(filter(lambda lf: d_val[corresp[lf]]>=0.01,rnlg))
    if (len(toz)>0): nd_val[toz]=0.
    # sum all contribution on owner
    nd_valv.scatter_reverse(la.InsertMode.add)
    # compute final value on owner into d_val (for owner dofs/mesh index are the same)
    if ( nloc>1):
        d_val[:nloc] = np.maximum(np.multiply(nd_val[:nloc] , num_edges_per_nodes[:nloc]), d_val[:nloc])
    # impose owner on ghost
    d.x.scatter_forward()
    #print(iter_smooth,"first creation",MPI.Wtime()-dstart)
    #dstart=MPI.Wtime()
    # prod again to smooth and enlarge
    nd_val[:]=mat.dot(d_val)
    # sum, compute, impose
    nd_valv.scatter_reverse(la.InsertMode.add)
    if ( nloc>1):
        d_val[:nloc] = np.maximum(np.multiply(nd_val[:nloc] , num_edges_per_nodes[:nloc]), d_val[:nloc])
    d.x.scatter_forward()
    #print(iter_smooth,"second creation",MPI.Wtime()-dstart)
dt[6]+=MPI.Wtime()-start # 4.2 Define damage

start=MPI.Wtime()
# specific field to output partition ========================
part=fem.Function(DGspace,name="partition")
part.x.array.fill(MPI.COMM_WORLD.Get_rank())
dt[12]+=MPI.Wtime()-start # 8 Outputs

# Elasticity ======================
start=MPI.Wtime()
# Elasticity parameters
# fill Young modulus field ==================================
libc = CDLL("libc.so.6") # Should be the same lib as the one used by the c++ programs to make 
                         # this python version the same as the C++ version
libc.srand(6575) # fix seed so that E_range is the same across process and for all runs

E_slop=(1.e8-5.e6)/199.
E_range=np.empty(200)
for i in range(200):
    E_range[i]=E_slop*(libc.rand()%200)+5.e+6
E=fem.Function(DGspace,name="E")
E.x.array[:]=E_range[cell_tag.values[:]%200]
nu=fem.Constant(domain,0.3)
# derived
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
dt[4]+=MPI.Wtime()-start # 4.1 Material constant
start=MPI.Wtime()
# functions =======================
du = TrialFunction(space) # Incremental displacement
delta_u = TestFunction(space) #Test function
u=fem.Function(space,name='disp_FEniCSx_symb_sym_py') # Displacement from previous iteration
# strain
def eps(v):
    return 0.5*(grad(v) + grad(v).T)
# symbolic ========================
# sympolic variables
l,m,dam,psi,al,al1,al2,e11,e22,e12=sp.symbols('lmbda mu d psi alpha alpha1 alpha2 eps11 eps22 eps12')
# put in place strain tensor as a symetric 2x2 matrix (eps21=eps12)
T = sp.Matrix([[e11,e12],[e12,e22]], symmetric=True, real=True)
# Compute the symbolic expression for eigenvalues by sympy
eigv = T.eigenvals(multiple=True) 
# potential
psi=(T.trace()**2)*(1-al*dam)*l/2+m*((1-al1*dam)*eigv[0]**2+(1-al2*dam)*eigv[1]**2)
# first derivative on strain => sigma
siga=sp.simplify(sp.diff(psi,T))
#sig=[siga[0,x] for x in [0,1]]
#sig.append(siga[1,1])
sig=[siga[0,x] for x in [0]]
sig.append(siga[0,1]/2.)
sig.append(siga[1,1])
# i1 strain invariant (trace)
i1=T.trace()
# i2 strain invariant (1/2(trace(eps^2)-trace(eps)^2)
i2=sp.simplify(((T*T).trace()-i1*i1))/2
# symbolic to string ==============
eig_expr = list(map(str, eigv))
sig_expr = list(map(str, sig))
i1_expr = str(i1)
i2_expr = str(i2)
#
# stress ==========================
def sigma_strain_not_nul(eps11,eps22,eps12,eps_i1):
    eig1,eig2=map(eval, eig_expr)
    alpha1 = conditional(lt(eig1,0.0),0.0,1.0)
    alpha2 = conditional(lt(eig2,0.0),0.0,1.0)
    alpha = conditional(lt(eps_i1,0.0),0.0,1.0)
    s00,s01,s11=map(eval, sig_expr)
    v= ([s00,s01],[s01,s11])
    return ufl.tensors.as_tensor(v)

def sig_sym_dam(strain,trace): 
    # TODO: in general case, for Jacobian, needs to limit d to avoid  matrix term 
    # that may lead to singular operator
    return (1.-d)*(2.0*mu*strain + lmbda*trace*Identity(2))

def sig_dam(T,eps11,eps22,eps12,eps_i1):
    eps_i2=eval(i2_expr)
    return  conditional(eq(eps_i1,0.0),conditional(eq(eps_i2,0.0),sig_sym_dam(T,eps_i1),sigma_strain_not_nul(eps11,eps22,eps12,eps_i1)),sigma_strain_not_nul(eps11,eps22,eps12,eps_i1))

def sigma(x): 
    T=eps(x)
    eps11=T[0,0]
    eps12=T[0,1]
    eps22=T[1,1]
    eps_i1=eval(i1_expr)
    return  conditional(gt(d,0.0),sig_dam(T,eps11,eps22,eps12,eps_i1),sig_sym_dam(T,eps_i1))
#
dt[9]+=MPI.Wtime()-start # 7.1 Nonlinear form creation
start=MPI.Wtime()
# Dirichelet ======================
bc_nodes=mesh.locate_entities(domain, dim0, lambda x: np.isclose(x[0], 0.))
bdofs=fem.locate_dofs_topological(space,dim0,bc_nodes)
bcl=fem.dirichletbc(np.array([0.,0.]),bdofs,space)
imp_nodes=mesh.locate_entities(domain, dim0, lambda x: np.isclose(x[0], 1.))
idofs=fem.locate_dofs_topological(space,dim0,imp_nodes)
if trac==True:
    imp=0.01
else:
    imp=-0.01
bcr=fem.dirichletbc(np.array([imp,0.]),idofs,space)
bcs=[bcl,bcr]
dt[7]+=MPI.Wtime()-start # 5.2 Dirichlet setting
# loads ===========================
start=MPI.Wtime()
# volumic force
f=fem.Function(space,name="f")
def f_vol_func(x):
    nb_row = len(x[1])
    assert nb_row == len(x[0])
    nb_col_v = 2
    v = np.empty((nb_col_v,nb_row))
    xf=100000.
    v[0,:]= -xf*np.power(x[0]-0.5,3)*(1600.*np.square(x[1]-0.5)-500.)
    v[1,:]= 0.
    return v
f.interpolate(f_vol_func)
#t = Constant(space)
#n = FacetNormal(domain)
dt[8] += MPI.Wtime() - start  # 5.2 Neuman setting
#
start=MPI.Wtime()
# quadrature ======================
metadata = {"quadrature_rule": "default", "quadrature_degree": 1}
dxx = measure.Measure("dx", domain = domain, metadata = metadata)
# Non linear function =============
F=inner(sigma(u), eps(delta_u))*dxx -  inner(f, delta_u)*dx  #- dot(t*n,delta_u)*ds(0)
# Jacobian of F ===================
J=ufl.derivative(F,u,du)
problem = fem.petsc.NonlinearProblem(F, u, bcs,J)
dt[9]+=MPI.Wtime()-start # 7.1 Nonlinear form creation
start=MPI.Wtime()
# Set non linear problem ==========
newton_solver = nls.petsc.NewtonSolver(MPI.COMM_WORLD, problem)
newton_solver.rtol = 1.e-7
newton_solver.atol = 5.e-8
ksp = newton_solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "cg"
opts[f"{option_prefix}ksp_rtol"] = "1e-12"
opts[f"{option_prefix}pc_type"] = "hypre"
opts[f"{option_prefix}pc_hypre_type"] = "boomeramg"
opts[f"{option_prefix}pc_hypre_type"] = "boomeramg"
opts[f"{option_prefix}pc_hypre_boomeramg_coarsen_type"] = "HMIS"
opts[f"{option_prefix}pc_hypre_boomeramg_relax_type_up"] = "l1scaled-SOR/Jacobi"
opts[f"{option_prefix}pc_hypre_boomeramg_relax_type_down"] = "l1scaled-SOR/Jacobi"
opts[f"{option_prefix}pc_hypre_boomeramg_relax_type_coarse"] = "Gaussian-elimination"
opts[f"{option_prefix}pc_hypre_boomeramg_interp_type"] = "ext+i"
opts[f"{option_prefix}pc_hypre_boomeramg_numfunctions"] = "2"
opts[f"{option_prefix}pc_hypre_boomeramg_agg_nl"] = "0"
opts[f"{option_prefix}pc_hypre_boomeramg_strong_threshold"] = "0.25"
#opts[f"{option_prefix}pc_hypre_boomeramg_print_debug"] = "1"
opts[f"{option_prefix}pc_hypre_boomeramg_print_statistics"] = "1"
ksp.setFromOptions()
dt[10]+=MPI.Wtime()-start # 7.2 Solver creation
start=MPI.Wtime()
# non linear resolution ===========
r=newton_solver.solve(u)
if (do_print): print(r)
dt[11]+=MPI.Wtime()-start # 7.3 NonLinear resolution
start=MPI.Wtime()
interpolation_point=DGVspace.element.interpolation_points()
# strain field ====================
def reps(v):
    re=eps(v)
    return ufl.as_vector([re[0,0],re[0,1],re[1,1]])

strain=fem.Function(DGVspace,name="strain_FeniCSx_symb_sym_py")
strain_exp=fem.Expression(reps(u), interpolation_point,MPI.COMM_WORLD)
#breakpoint()
strain.interpolate(strain_exp)

# stress field ====================
def rsig(v):
    re=sigma(v)
    return ufl.as_vector([re[0,0],re[0,1],re[1,1]])

stress=fem.Function(DGVspace,name="stress_FeniCSx_symb_sym_py")
stress_exp=fem.Expression(rsig(u), interpolation_point,MPI.COMM_WORLD)
#breakpoint()
stress.interpolate(stress_exp)
dt[15]+=MPI.Wtime()-start # 8.1 strain/stress computation

# outputs =========================
start=MPI.Wtime()
filename = Path(f"run/asym_elasto_damage_trac_vol_symb_sym_scal_py_{MPI.COMM_WORLD.Get_size():d}.bp")
with VTXWriter(MPI.COMM_WORLD, filename, [d]) as ofile:
    ofile.write(0.)
filename = Path(f"run/asym_elasto_damage_trac_vol_symb_sym_DGscal_py_{MPI.COMM_WORLD.Get_size():d}.bp")
with VTXWriter(MPI.COMM_WORLD, filename, [E,part]) as ofile:
    ofile.write(0.)
if trac:
    filename = Path(f"run/asym_elasto_damage_trac_vol_symb_sym_tens_py_{MPI.COMM_WORLD.Get_size():d}.bp")
else:
    filename = Path(f"run/asym_elasto_damage_comp_vol_symb_sym_tens_py_{MPI.COMM_WORLD.Get_size():d}.bp")
with VTXWriter(MPI.COMM_WORLD, filename, [strain,stress]) as ofile:
    ofile.write(0.)
if trac:
    filename = Path(f"run/asym_elasto_damage_trac_vol_symb_sym_vect_py_{MPI.COMM_WORLD.Get_size():d}.bp")
else:
    filename = Path(f"run/asym_elasto_damage_comp_vol_symb_sym_vect_py_{MPI.COMM_WORLD.Get_size():d}.bp")
with VTXWriter(MPI.COMM_WORLD, filename, [u,f]) as ofile:
    ofile.write(0.)
dt[12]+=MPI.Wtime()-start # 8 Outputs
# measure output ==================
dt[0]=MPI.Wtime()-start_all # All
def fmt_out(i,m):
    print(f"| {min_dt[i]:12.5f} | {max_dt[i]:12.5f} | {dt[i]:12.5f} | {(100 * dt[i] / dt[0]):12.1f} | {(avg_dt[i] / nb_proc):12.5f} | {(100 * avg_dt[i] / avg_dt[0]):12.1f} | {m:>42} |")
if  MPI.COMM_WORLD.Get_rank()>0:
    max_dt=np.zeros(1)
    min_dt=np.zeros(1)
    avg_dt=np.zeros(1)
else:
    max_dt=np.zeros(16)
    min_dt=np.zeros(16)
    avg_dt=np.zeros(16)
MPI.COMM_WORLD.Reduce(dt, max_dt, op=MPI.MAX, root=0)
MPI.COMM_WORLD.Reduce(dt, min_dt, op=MPI.MIN, root=0)
MPI.COMM_WORLD.Reduce(dt, avg_dt, op=MPI.SUM, root=0)
#print(dt)
if do_print:
    print(f"{'='*136}")
    print(f"| {'val min':>12} | {'val max':>12} | {'val':>12} | {'% total':>12} | {'val avg':>12} | {'% total avg':>12} | {' ':>42} |")
    fmt_out(0, "All")
    fmt_out(1, "1. Initialize ")
    fmt_out(2, "2.1 Read the mesh")
    fmt_out(3, "2.2 Refine the mesh")
    fmt_out(5, "3.1  Define space")
    fmt_out(6, "3.2  Define damage")
    fmt_out(4, "4.1  Material constant")
    fmt_out(7, "5.1 Dirichlet setting")
    fmt_out(8, "5.2 Neuman setting")
    #fmt_out(13, "6.3 Create and assemble elementary vector")
    #fmt_out(14, "6.4 Create and assemble elementary matrix")
    fmt_out(9, "7.1 Nonlinear form creation")
    fmt_out(10, "7.2 Solver creation")
    fmt_out(11, "7.3 NonLinear resolution")
    fmt_out(12, "8 Outputs")
    fmt_out(15, "8.1 strain/stress computation")
    print(f"{'='*136}")
