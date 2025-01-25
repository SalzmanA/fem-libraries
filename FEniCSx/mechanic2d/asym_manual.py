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
# stress computation ==============
def sig_sym_dam(strain,trace): 
    return (1.-d)*(2.0*mu*strain + lmbda*trace*Identity(2))

def sig_dam_prod(vect,diag): 
    ed=diag[0]*vect[0][0]*vect[1][0]+diag[1]*vect[1][1]*vect[0][1]
    v=([diag[0]*vect[0][0]*vect[0][0]+diag[1]*vect[0][1]*vect[0][1],ed],[ed,diag[0]*vect[1][0]*vect[1][0]+diag[1]*vect[1][1]*vect[1][1]])
    return ufl.tensors.as_tensor(v)

def sig_dam_eignv1(sl): 
    evecl= [[1.,0.],[0.,1.]]
    return sig_dam_prod(evecl,sl)

def sig_dam_eignv2(strai,e1,e2,sl): 
    v1=e1-strai[1,1]
    v2=e2-strai[1,1]
    v3=strai[1,0]
    v32=v3*v3
    n1=v1*v1+v32
    n1=sqrt(n1)
    n2=v2*v2+v32
    n2=sqrt(n2)
    evecl= [[v1/n1,v2/n2],[v3/n1,v3/n2]]
    return sig_dam_prod(evecl,sl)

def sig_dam_eignv(stra,ei1,ei2,a1,a2,a): 
    factor=2.*mu
    gamma=0.5*lmbda/mu
    c=1.-a*d
    c1=1.-a1*d
    c2=1.-a2*d
    D0=factor*(c1+gamma*c)
    D1=factor*gamma*c
    D2=factor*(c2+gamma*c)
    evall=[D0*ei1+D1*ei2,D1*ei1+D2*ei2]
    return conditional(lt(stra[1,0],1.e-12),conditional(gt(stra[1,0],-1.e-12),sig_dam_eignv1(evall),sig_dam_eignv2(stra,ei1,ei2,evall)),sig_dam_eignv2(stra,ei1,ei2,evall))

def sig_dam_strain(straint,iv1,iv2): 
    delta=iv1*iv1+4*iv2
    r=sqrt(delta)
    eig1 = (iv1 + r) / 2.
    eig2 = (iv1 - r) / 2.
    alpha1 = conditional(lt(eig1,0.0),0.0,1.0)
    alpha2 = conditional(lt(eig2,0.0),0.0,1.0)
    alpha = conditional(lt(eig1+eig2,0.0),0.0,1.0)
    return  conditional(eq(d,1.0), conditional(eq(alpha1,1.0), conditional(eq(alpha2,1.0), conditional(eq(alpha,1.0), sig_sym_dam(straint,iv1),sig_dam_eignv(straint,eig1,eig2,alpha1,alpha2,alpha)), sig_dam_eignv(straint,eig1,eig2,alpha1,alpha2,alpha)), sig_dam_eignv(straint,eig1,eig2,alpha1,alpha2,alpha)), sig_dam_eignv(straint,eig1,eig2,alpha1,alpha2,alpha))

def sig_dam(strain,inv1): 
    i2=0.5*(tr(strain*strain)-inv1*inv1)
    return  conditional(eq(inv1,0.0),conditional(eq(i2,0.0),sig_sym_dam(strain,inv1),sig_dam_strain(strain,inv1,i2)),sig_dam_strain(strain,inv1,i2))

def sigma(x): 
    T=eps(x);
    i1=tr(T)
    return  conditional(gt(d,0.0),sig_dam(T,i1),sig_sym_dam(T,i1))
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
