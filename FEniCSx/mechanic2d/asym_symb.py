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
# symbolic ========================
# sympolic variables
l,m,dam,psi,al,al1,al2=sp.symbols('lmbda mu d psi alpha alpha1 alpha2')
# put in place strain tensor as a symetric 2x2 matrix (symetric=True argument) 
T = sp.Matrix(2, 2, lambda i, j: sp.Symbol('T[%d, %d]' % (i, j), symmetric=True, real=True))
# Compute the symbolic expression for eigenvalues by sympy
eigv = T.eigenvals(multiple=True) 
# potential
psi=(T.trace()**2)*(1-al*dam)*l/2+m*((1-al1*dam)*eigv[0]**2+(1-al2*dam)*eigv[1]**2)
# first derivative on strain => sigma
sig=sp.simplify(sp.diff(psi,T))
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
def sigma_strain_not_nul(T,tre):
    eig1,eig2=map(eval, eig_expr)
    alpha1 = conditional(lt(eig1,0.0),0.0,1.0)
    alpha2 = conditional(lt(eig2,0.0),0.0,1.0)
    alpha = conditional(lt(tre,0.0),0.0,1.0)
    s00,s01,s10,s11=map(eval, sig_expr)
    v= ([s00,s01],[s10,s11])
    return ufl.tensors.as_tensor(v)

def sig_sym_dam(strain,trace): 
    # TODO: in general case, for Jacobian, needs to limit d to avoid  matrix term 
    # that may lead to singular operator
    return (1.-d)*(2.0*mu*strain + lmbda*trace*Identity(2))

def sig_dam(T,eps_i1):
    eps_i2=eval(i2_expr)
    return  conditional(eq(eps_i1,0.0),conditional(eq(eps_i2,0.0),sig_sym_dam(T,eps_i1),sigma_strain_not_nul(T,eps_i1)),sigma_strain_not_nul(T,eps_i1))

def sigma(x): 
    T=eps(x);
    eps_i1=eval(i1_expr)
    return  conditional(gt(d,0.0),sig_dam(T,eps_i1),sig_sym_dam(T,eps_i1))
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
