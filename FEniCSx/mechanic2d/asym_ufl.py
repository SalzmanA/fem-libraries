import sympy as sp
import ufl
from ufl.classes import (Mesh,FunctionSpace,Coefficient,Constant,
                         TrialFunction,TestFunction,FacetNormal)
from ufl import ( sqrt,grad,inner,tr,Identity,dot,dx,ds,
                 conditional,lt,gt,eq,measure)
import basix


# element & space =================
elem = basix.ufl.element("Lagrange", "triangle", 1, shape=(2, ))
element_scal = basix.ufl.element("Lagrange", "triangle", 1)
DGelement = basix.ufl.element("DG", "triangle", 0)
domain = Mesh(elem)
space = FunctionSpace(domain, elem)
DGspace = FunctionSpace(domain, DGelement)
space_scal = FunctionSpace(domain, element_scal)

# damage ==========================
d=Coefficient(space_scal)
# Elasticity ======================
# Elasticity parameters
E=Coefficient(DGspace)
nu=Constant(space)
# derived
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
# functions =======================
du = TrialFunction(space) # Incremental displacement
delta_u = TestFunction(space) #Test function
u = Coefficient(space)  # Displacement from previous iteration
# strain ==========================
def eps(v):
    return 0.5*(grad(v) + grad(v).T)
# potential =======================

def psi_sym_dam(T,i1):
    return (1 - d) * (i1 * i1 * lmbda / 2. + mu * (T[0,0] * T[0,0] + T[1,1] * T[1,1] + T[1,0] * T[1,0] + T[0,1] * T[0,1]))

def psi_non_null_strain(T,i1,i2):
    delta=i1*i1+4*i2
    r=sqrt(delta)
    eig1 = (i1 + r) / 2.
    eig2 = (i1 - r) / 2.
    alpha1 = conditional(lt(eig1,0.0),0.0,1.0)
    alpha2 = conditional(lt(eig2,0.0),0.0,1.0)
    alpha = conditional(lt(eig1+eig2,0.0),0.0,1.0)
    return i1*i1*(1. - alpha * d) * lmbda / 2. + mu * ((1 - alpha1 * d) * eig1 * eig1 + (1. - alpha2 * d) * eig2 * eig2)

def psi(S):
    i1=tr(S)
    i2=0.5*(tr(S*S)-i1*i1)
    return  conditional(eq(i1,0.0),conditional(eq(i2,0.0),psi_sym_dam(S,i1),psi_non_null_strain(S,i1,i2)),psi_non_null_strain(S,i1,i2))
#
# stress ==========================
def sig_sym_dam(T): 
    T=ufl.variable(T)
    i1=tr(T)
    v= ufl.diff(psi_sym_dam(T,i1),T)
    return ufl.tensors.as_tensor(v)

def sig_dam(T): 
    T=ufl.variable(T)
    v= ufl.diff(psi(T),T)
    return ufl.tensors.as_tensor(v)

def sigma(x): 
    T=eps(x);
    return  conditional(gt(d,0.0),sig_dam(T),sig_sym_dam(T))
#
# loads ===========================
# volumic force
f = Coefficient(space)
# Traction force
#t = Constant(space)
#n = FacetNormal(domain)
# quadrature ======================
metadata = {"quadrature_rule": "default", "quadrature_degree": 1}
dxx = measure.Measure("dx", domain = domain, metadata = metadata)
# Non linear function =============
F=inner(sigma(u), eps(delta_u))*dxx -  inner(f, delta_u)*dx  #- dot(t*n,delta_u)*ds(0)
# Jacobian of F ===================
J=ufl.derivative(F,u,du)
# Expression ======================
# reduced strain expression
def reps(v):
    re=eps(v)
    return ufl.as_vector([re[0,0],re[0,1],re[1,1]])
def rsig(v):
    re=sigma(v)
    return ufl.as_vector([re[0,0],re[0,1],re[1,1]])
def energ(v):
    strain=eps(v)
    stress=sigma(v)
    return inner(strain,stress)
center=1./3.
expressions=[(reps(u), [[center, center]]),(rsig(u), [[center, center]]),(energ(u), [[center, center]])]
