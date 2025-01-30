/// \file
#include <fstream>
#include <iomanip>
#include <iostream>

#define MAX_DAM 1.
//#define MAX_DAM 0.
//#define PETSC_SOLV
#define PETSC_LIKE_BOOMERANG_SETTING
//#define TEST_ADAPT_LINRTOL
//#define DEBUG_SQUARE
#define USE_TRAC
#define USE_VOLUME
#define USE_B
//#define USE_F_OUTSIDE
//#define USE_AD
#define USE_ADIOS_FOR_OUTPUT
//#define OUT_COMP

#include "mfem.hpp"
#ifdef USE_AD
#include "autodiff/admfem.hpp"
#endif

#define USE_ECST

//#define MFEM_THREAD_SAFE

#ifdef PETSC_SOLV
#ifndef MFEM_USE_PETSC
#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
#endif
#error Actual implementation is not working with PETSc as linear solver
#endif
#define SL(i, m)                                                                                                              \
   cout << "| " << fixed << setprecision(5) << setw(12) << min_measure[i] << " | " << setw(12) << max_measure[i] << " | "     \
        << setw(12) << sqrt(ecart_type_measure[i] * nb_proc - avg_measure[i] * avg_measure[i]) / nb_proc << " | " << setw(12) \
        << 100 * sqrt(ecart_type_measure[i] * nb_proc - avg_measure[i] * avg_measure[i]) / avg_measure[0] << " | "            \
        << setprecision(5) << setw(12) << avg_measure[i] / nb_proc << " | " << setprecision(1) << setw(12)                    \
        << 100 * avg_measure[i] / avg_measure[0];                                                                             \
   cout << " | " << setw(42) << m << " |" << endl;

#define NB_MEASURE 16
#ifdef USE_TRAC
#ifdef USE_ECST
#undef USE_ECST
#endif
#endif

#ifdef OUT_COMP
#ifdef USE_AD
#define IN_COMP
#undef OUT_COMP
#endif
#endif

#ifdef OUT_COMP
#ifdef USE_F_OUTSIDE
#define IN_COMP
#undef OUT_COMP
#endif
#endif

using namespace std;
using namespace mfem;

#ifdef IN_COMP
#include <algorithm>
#include <fstream>
#endif

#ifdef USE_F_OUTSIDE
// simple derivation of original VectorDomainLFIntegrator to add profiling
class VectorDomainLFIntegratorProf : public VectorDomainLFIntegrator
{
  private:
   double measure;

  public:
   VectorDomainLFIntegratorProf(VectorCoefficient &QF) : VectorDomainLFIntegrator(QF), measure(0.) {}

   virtual void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
   {
      const double start = MPI_Wtime();
      VectorDomainLFIntegrator::AssembleRHSElementVect(el, Tr, elvect);
      measure += MPI_Wtime() - start;
   }
   double getMeasure() { return measure; }
};
#endif

#ifdef USE_AD
/// Functor used in automatic differentiation that compute the potential:
// psi=(strain.trace()**2)*(1-alpha*dam)*l/2+m*((1-alpha1*dam)*eigv[0]**2+(1-alpha2*dam)*eigv[1]**2)
// (i.e. asymmetric traction/compression damaged elasticity potential)
/// State vector strain store strain tensor in the following order:
///  0     1      2      3
/// eps11, eps21, eps12, eps22
/// and parameters provided in vector vparam are Lame and damage
template <typename TDataType, typename TParamVector, typename TStateVector, int state_size, int param_size>
class Potential
{
   double limit = 1.e-12;
   double mlimit = -1.e-12;

  public:
   TDataType operator()(TParamVector &vparam, TStateVector &strain)
   {
      MFEM_ASSERT(state_size == 4, "Potential state_size should be equal to 4!");
      MFEM_ASSERT(param_size == 3, "Potential param_size should be equal to 3!");
      real_t d = vparam[2];  // Damage

      // computing the invariants
      TDataType I1 = strain[0] + strain[3];
      TDataType I2 = strain[1] * strain[2] - strain[0] * strain[3];

      // non null tensor
      if (I1 > limit || I2 > limit || I1 < mlimit || I2 < mlimit)
      {
         // computing the eigen values
         TDataType delta = I1 * I1 + 4. * I2;
         MFEM_ASSERT(delta > -1.e-12, "!!! Strain tensor not correctly defined: rank deficient !!!")
         TDataType r = sqrt(delta);
         TDataType ev1 = (I1 + r) / 2.;
         TDataType ev2 = (I1 - r) / 2.;

         // alpha's
         real_t alpha, alpha1, alpha2;
         if (ev1 >= 0)
            alpha1 = 1.;
         else
            alpha1 = 0.;
         if (ev2 >= 0)
            alpha2 = 1.;
         else
            alpha2 = 0.;
         if ((ev1 + ev2) >= 0)
            alpha = 1.;
         else
            alpha = 0.;

         // return the asymetric potential
         return I1 * I1 * (1. - alpha * d) * vparam[0] / 2. +
                vparam[1] * ((1 - alpha1 * d) * ev1 * ev1 + (1. - alpha2 * d) * ev2 * ev2);
      }
      // null tensor
      else
      {
         // switch to simple linear potential to avoid singular second derivative (due to sqr at zero). Even if potential is
         // defacto null in this case it's derivatives are not. Full expretion so that AD works well.
         return (1 - d) * (I1 * I1 * vparam[0] / 2. + vparam[1] * (strain[0] * strain[0] + strain[3] * strain[3] +
                                                                   strain[1] * strain[1] + strain[2] * strain[2]));
      }
   }
};

// function used to compute asymetric stress tensor
void asym_stress(const FiniteElement &el, ElementTransformation &Tr, const IntegrationPoint &ip, const real_t &w,
                 const DenseMatrix &disp, QuadratureFunctionCoefficient &dam, Coefficient &lambda, Coefficient &mu, DenseMatrix &dshape,
                 DenseMatrix &gdshape, DenseMatrix &strain, Vector &vparam,
                 mfem::QFunctionAutoDiff<Potential, 4, 3> &ad_potentiel, DenseMatrix &sig)
{
   // get damage
   real_t d = dam.Eval(Tr, ip);
   // d=0.;
   // get lambda
   real_t l = lambda.Eval(Tr, ip);
   // get mu
   real_t m = mu.Eval(Tr, ip);

   // get gradiant of form function
   el.CalcDShape(ip, dshape);
   // pass it in physical space
   Mult(dshape, Tr.InverseJacobian(), gdshape);
   // computes gradiant of displacement
   MultAtB(disp, gdshape, strain);
   // set as strain
   strain.Symmetrize();

   // if non null damage
   if (d > 0.)
   {
      vparam[0] = l * w;
      vparam[1] = m * w;
      vparam[2] = d;
      Vector vstrain(strain.GetData(), 4);
      Vector vsig(sig.GetData(), 4);
      ad_potentiel.Grad(vparam, vstrain, vsig);
   }
   // if null damage
   else
   {
      // simple linear stress: integration of E:strain
      // multiply here by weight
      real_t m2plw = w * (2 * m + l);
      real_t lw = l * w;
      sig(0, 0) = m2plw * strain(0, 0) + lw * strain(1, 1);
      sig(1, 1) = m2plw * strain(1, 1) + lw * strain(0, 0);
      sig(1, 0) = sig(0, 1) = w * m * (strain(0, 1) + strain(1, 0));
   }

   // cout << "stress" << endl;
   // sig.Print(cout);
}
#else  // non AD implementation
// function used to compute asymetric stress tensor
void asym_stress(const FiniteElement &el, ElementTransformation &Tr, const IntegrationPoint &ip, const real_t &w,
                 const DenseMatrix &disp, QuadratureFunctionCoefficient &dam, Coefficient &lambda, Coefficient &mu, DenseMatrix &dshape,
                 DenseMatrix &gdshape, DenseMatrix &strain, real_t &limit, real_t &mlimit, Vector &eval, Vector &eigns,
                 DenseMatrix &evec, DenseMatrix &sig)
{
   // get damage
   real_t d = dam.Eval(Tr, ip);
   // d=0.;
   // get lambda
   real_t l = lambda.Eval(Tr, ip);
   // get mu
   real_t m = mu.Eval(Tr, ip);

   // get gradiant of form function
   el.CalcDShape(ip, dshape);
   // pass it in physical space
   Mult(dshape, Tr.InverseJacobian(), gdshape);
   // computes gradiant of displacement
   MultAtB(disp, gdshape, strain);
   // set as strain
   strain.Symmetrize();

   // if non null damage
   if (d > 0.)
   {
      // computing the invariants
      real_t I1 = strain(0, 0) + strain(1, 1);
      real_t I2 = strain(0, 1) * strain(0, 1) - strain(0, 0) * strain(1, 1);

      // if non null tensor
      if (I1 > limit || I2 > limit || I1 < mlimit || I2 < mlimit)
      {
         // computing the eigen values
         real_t delta = I1 * I1 + 4 * I2;
         MFEM_ASSERT(delta > mlimit, "!!! Strain tensor not correctly defined: rank deficient !!!")
         real_t r = sqrt(max(0., delta));
         // cout << "r " << r << " delta " << delta << endl;
         eval[0] = (I1 + r) / 2.;
         eval[1] = (I1 - r) / 2.;

         // cout << "eign..." << endl;

         // modifying values of alpha's based on the signs of eigenvalues
         real_t alpha, alpha1, alpha2;
         if (eval[0] >= 0)
            alpha1 = 1.;
         else
            alpha1 = 0.;
         if (eval[1] >= 0)
            alpha2 = 1.;
         else
            alpha2 = 0.;
         if ((eval[0] + eval[1]) >= 0)
            alpha = 1.;
         else
            alpha = 0.;

         // genral case : not pure traction and full damage (i.e. alpha_i=1=damage)
         if (!((d == 1.) && (alpha == 1) && (alpha1 == 1) && (alpha2 == 1)))
         {
            // computing the eigenvectors
            if (fabs(strain(1, 0)) > limit)
            {
               evec(0, 0) = eval[0] - strain(1, 1);
               evec(0, 1) = eval[1] - strain(1, 1);
               evec(1, 0) = evec(1, 1) = strain(1, 0);
               real_t norms[2];
               evec.Norm2(norms);
               evec(0, 0) /= norms[0];
               evec(1, 0) /= norms[0];
               evec(0, 1) /= norms[1];
               evec(1, 1) /= norms[1];
            }
            else
            {
               evec(0, 0) = evec(1, 1) = 1.;
               evec(1, 0) = evec(0, 1) = 0.;
            }
            real_t temp = 2. * m * w;    //  E / (1 + nu) multiply here by weight
            real_t gamma = 0.5 * l / m;  // nu / (1 - 2 * nu)

            real_t c = 1 - alpha * d;
            real_t c1 = 1 - alpha1 * d;
            real_t c2 = 1 - alpha2 * d;
            real_t D0 = temp * (c1 + gamma * c);
            real_t D1 = temp * gamma * c;
            real_t D2 = temp * (c2 + gamma * c);
            eigns[0] = D0 * eval[0] + D1 * eval[1];
            eigns[1] = D1 * eval[0] + D2 * eval[1];
            // cout << "eigns" << endl;
            // eigns.Print(cout);
            // cout << "evec" << endl;
            // evec.Print(cout);

            MultADAt(evec, eigns, sig);
         }
         // pure traction and full damage: alpha_i=1=damage
         else
         {
            sig = 0.;  // stress=0 as
         }
      }
      // null strain tensor
      else
      {
         sig = 0.;  // stress=0
      }
   }
   // if null damage
   else
   {
      // simple linear stress: integration of E:strain
      // multiply here by weight
      real_t m2plw = w * (2 * m + l);
      real_t lw = l * w;
      sig(0, 0) = m2plw * strain(0, 0) + lw * strain(1, 1);
      sig(1, 1) = m2plw * strain(1, 1) + lw * strain(0, 0);
      sig(1, 0) = sig(0, 1) = w * m * (strain(0, 1) + strain(1, 0));
   }

   // cout << "stress" << endl;
   // sig.Print(cout);
}
#endif

// to project strain on vectorial DGspace for output
class strainTensor : public VectorCoefficient
{
  public:
   strainTensor(ParGridFunction &u_) : VectorCoefficient(3), u(u_), nd(3), dim(2) {}
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
#ifdef MFEM_THREAD_SAFE
      DenseMatrix eps;
#endif
      // strain computation
      u.GetVectorGradient(T, eps);
      eps.Symmetrize();
      // store in V:
      // eps11 in x, eps12 in y, eps12=eps22 in z
      V.Elem(0) = eps(0, 0);
      V.Elem(1) = eps(0, 1);
      V.Elem(2) = eps(1, 1);
   }
   void Eval(DenseMatrix &M, ElementTransformation &T, const IntegrationRule &ir) override
   {
      cout << "Not implemented !!" << endl;
      throw -345;
   }

  private:
   ParGridFunction &u;
   const int nd, dim;
#ifndef MFEM_THREAD_SAFE
   DenseMatrix eps;
#endif
};

// to project stress on vectorial DGspace for output
class stressTensor : public VectorCoefficient
{
  public:
   stressTensor(ParGridFunction &u_, QuadratureFunctionCoefficient &dam_, Coefficient &l, Coefficient &m)
       : VectorCoefficient(3),
         u(u_),
         nd(3),
         dim(2),
         dam(dam_),
         lambda(l),
         mu(m),
         one(1.)
#ifndef USE_AD
         ,
         limit(1.e-12),
         mlimit(-1.e-12)
#endif
#ifndef MFEM_THREAD_SAFE
         ,
         sig(dim),
         dshape(nd, dim),
         gdshape(nd, dim),
         strain(dim),
#ifdef USE_AD
         vparam(3)
#else
         eval(dim),
         eigns(dim),
         evec(dim)
#endif
#endif
   {
   }
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
#ifdef MFEM_THREAD_SAFE
      DenseMatrix sig(dim), dshape(nd, dim), gdshape(nd, dim), strain(dim), disp;
      Vector lval;
#endif
      // get element (assertion: same spacefor u,dam,lambda and mu
      const FiniteElement *FElem = u.ParFESpace()->GetFE(T.ElementNo);
      // get disp
      u.GetElementDofValues(T.ElementNo, lval);
      disp.UseExternalData(lval.GetData(), nd, dim);
      // stress computation
#ifdef USE_AD
#ifdef MFEM_THREAD_SAFE
      Vector vparam(3);
#endif
      asym_stress(*FElem, T, ip, one, disp, dam, lambda, mu, dshape, gdshape, strain, vparam, ad_potentiel, sig);
#else
#ifdef MFEM_THREAD_SAFE
      DenseMatrix evec(dim);
      Vector eval(dim), eigns(dim);
#endif
      asym_stress(*FElem, T, ip, one, disp, dam, lambda, mu, dshape, gdshape, strain, limit, mlimit, eval, eigns, evec, sig);
#endif
      // store in V:
      // sig11 in x, sig12 in y, sig12=sig22 in z
      V.Elem(0) = sig(0, 0);
      V.Elem(1) = sig(0, 1);
      V.Elem(2) = sig(1, 1);
   }
   void Eval(DenseMatrix &M, ElementTransformation &T, const IntegrationRule &ir) override
   {
      cout << "Not implemented !!" << endl;
      throw -345;
   }

  private:
   ParGridFunction &u;
   const int nd, dim;
   QuadratureFunctionCoefficient &dam;
   Coefficient &lambda;
   Coefficient &mu;
   real_t one;
#ifdef USE_AD
   mfem::QFunctionAutoDiff<Potential, 4, 3> ad_potentiel;
#else
   real_t limit, mlimit;
#endif
#ifndef MFEM_THREAD_SAFE
   DenseMatrix sig, dshape, gdshape, strain, disp;
   Vector lval;
#ifdef USE_AD
   Vector vparam;
#else
   Vector eval, eigns;
   DenseMatrix evec;
#endif
#endif
};

#ifdef IN_COMP
class energyError : public Coefficient
{
  public:
   energyError(strainTensor &strain_, stressTensor &stress_) : strain(strain_), stress(stress_), eps(3), sig(3), sum(0.) {}
   real_t getSum() { return sum; }
   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
#ifdef MFEM_THREAD_SAFE
#error "not thread safe due to sum"
#endif
      // strain
      strain.Eval(eps, T, ip);
      // strain
      stress.Eval(sig, T, ip);
      // error
      real_t err = eps.Elem(0) * sig.Elem(0) + eps.Elem(1) * sig.Elem(1) + 2. * eps.Elem(2) * sig.Elem(2);
      sum += err;
      return err;
   }

  private:
   strainTensor &strain;
   stressTensor &stress;
   Vector eps, sig;
   real_t sum;
};
#endif

// Non linear integrator to compute an asymmetric traction/compression damaged elasticity law:
// 2D plane strain implementation
class damIntegrator : public NonlinearFormIntegrator
{
  public:
   damIntegrator(Coefficient &l, Coefficient &m, QuadratureFunctionCoefficient &d, const IntegrationPoint &ip_,
                 const IntegrationRule *ir
#ifdef USE_F_OUTSIDE
                 )
#else
                 ,
                 VectorQuadratureFunctionCoefficient &load_)
#endif
       : NonlinearFormIntegrator(ir),
         lambda(l),
         mu(m),
         dam(d),
         ip(ip_),
#ifndef USE_F_OUTSIDE
         load(load_),
#endif
         time_vect(0.),
         time_grad(0.),
         nd(3),
         dim(2),
         limit(1.e-12),
         mlimit(-1.e-12)
#ifndef MFEM_THREAD_SAFE
         ,
         dshape(nd, dim),
         strain(dim),
         sig(dim),
         hook(3),
         block(nd),
#ifdef USE_B
         B(nd * dim, 3),
         C(nd * dim, 3),
#else
         gxgxt(nd),
         gygyt(nd),
         gxgyt(nd),
         gygxt(nd),
         gxygxyt(nd),
#endif
#ifdef USE_AD
         hess(4, 4),
         vparam(3),
#else
         evec(dim),
         potentiel_sderivative(3),
         dedeps(3),
         M(3),
         eval(dim),
         eigns(dim),
#endif
         vl(dim * nd),
         shape(nd),
         f(dim)
#endif
   {
#ifndef MFEM_THREAD_SAFE
      gdshape.UseExternalData(vl.GetData(), nd, dim);
#endif
   }
   real_t GetElementEnergy(const FiniteElement &el, ElementTransformation &Ttr, const Vector &elfun)
   {
      cout << "in GetElementEnergy" << endl;
      throw -1.;  // to code
      return 0.;
   }

   void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
   {
      // cout << "in AssembleElementVector ===================================================" << endl;
      MFEM_ASSERT(nd == el.GetDof(), "damIntegrator tuned to work with triangle first order element")
      MFEM_ASSERT(dim == el.GetDim(), "damIntegrator tuned to work with triangle in 2D spatial space")
      MFEM_ASSERT(dim == Tr.GetSpaceDim(), "damIntegrator tuned to work with 2D vectorial field")
      double start = MPI_Wtime();

      // reshape vector result
      elvect.SetSize(nd * dim);

      // local matrices and vector
#ifdef MFEM_THREAD_SAFE
      // local allocation if thread safety is required
      DenseMatrix dshape(nd, dim), gdshape(nd, dim), strain(dim), sig(dim), disp, res;
      Vector shape(nd), f(dim);
#ifdef USE_AD
      Vector vparam(3);
      mfem::QFunctionAutoDiff<Potential, 4, 3> ad_potentiel;
#else
      DenseMatrix evec(dim);
      Vector eval(dim);
      Vector eigns(dim);
#endif
#endif

      // set vectors into dense matrix used as warper
      disp.UseExternalData(elfun.GetData(), nd, dim);
      res.UseExternalData(elvect.GetData(), nd, dim);

      elvect = 0.;
      // get weight
      Tr.SetIntPoint(&ip);
      const real_t wt = Tr.Weight();
      real_t w = ip.weight * wt;
#ifdef USE_AD
      asym_stress(el, Tr, ip, w, disp, dam, lambda, mu, dshape, gdshape, strain, vparam, ad_potentiel, sig);
#else
      asym_stress(el, Tr, ip, w, disp, dam, lambda, mu, dshape, gdshape, strain, limit, mlimit, eval, eigns, evec, sig);
#endif

      // compute sig:grad
      AddMult(gdshape, sig, res);
      /*if ((i == 0 && nip == 1) || (i == 2 && nip > 1))
      {
         cout << "res " << endl;
         res.Print(cout);
      }*/

#ifdef USE_VOLUME
#ifndef USE_F_OUTSIDE
      // Load contribution
      //===================
      // integration loop
      const int nip = IntRule->GetNPoints();
      for (int i = 0; i < nip; i++)
      {
         const IntegrationPoint &ipl = IntRule->IntPoint(i);
         Tr.SetIntPoint(&ipl);
         // Tr.Weight constant for order 1 element
         //real_t wl = ipl.weight * Tr.Weight();
         real_t wl = ipl.weight * wt;
         
         // get volumic load
         load.Eval(f, Tr, ipl);

         // cout << "f" << endl;
         // f.Print(cout);

         // get form function in physical space
         el.CalcPhysShape(Tr, shape);
         // add -f:v
         AddMult_a_VWt(-wl, shape, f, res);
      }
#endif
#endif

      time_vect += MPI_Wtime() - start;
   }

   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
   {
      // cout << "in AssembleElementGrad =====================================================" << endl;
      MFEM_ASSERT(nd == el.GetDof(), "damIntegrator tuned to work with triangle first order element")
      MFEM_ASSERT(dim == el.GetDim(), "damIntegrator tuned to work with 2D vectorial field")
      double start = MPI_Wtime();

      // reset result matrix
      elmat.SetSize(nd * dim);
      elmat = 0.0;

      // local matrices and vector
#ifdef MFEM_THREAD_SAFE
      // local allocation if thread safety is required
      DenseMatrix dshape(nd, dim), gdshape, strain(dim), hook(3, 3), block(nd), disp;
      Vector vl(dim * nd);
      gdshape.UseExternalData(vl.GetData(), nd, dim);
#ifdef USE_B
      DenseMatrix B(nd * dim, 3), C(nd * dim, 3), mgdshapex, mgdshapey;
      B = 0.;
#else
      DenseMatrix gxgxt(nd), gygyt(nd), gxgyt(nd), gygxt(nd), gxygxyt(nd);
#endif
#ifdef USE_AD
      DenseMatrix hess(4, 4);
      Vector vparam(3);
      mfem::QFunctionAutoDiff<Potential, 4, 3> ad_potentiel;
#else
      DenseMatrix potentiel_sderivative(3), dedeps(3), M(3);
#endif
#endif
      // warping of vl from local or member depending on thread safety
      Vector gdshapex(vl.GetData(), nd);
      Vector gdshapey(vl.GetData() + nd, nd);
      disp.UseExternalData(elfun.GetData(), nd, dim);
#ifdef USE_B
      // warping of gdshapex,gdshapey local but warper is local or member depending on thread safety
      mgdshapex.UseExternalData(gdshapex.GetData(), nd, 1);
      mgdshapey.UseExternalData(gdshapey.GetData(), nd, 1);
#endif
#ifdef USE_AD
      Vector vstrain(strain.GetData(), dim * dim);
#endif

      // get weight
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      // get lambda
      real_t l = lambda.Eval(Tr, ip);
      // get mu
      real_t m = mu.Eval(Tr, ip);

      // get gradiant of form function
      el.CalcDShape(ip, dshape);
      // cout << "dshape" << endl;
      // dshape.Print(cout);
      // pass it in physical space
      Mult(dshape, Tr.InverseJacobian(), gdshape);
      // cout << "gdshape" << endl;
      // gdshape.Print(cout);
#ifdef USE_B
      // set B (always same block thus zero block remain unchanged thus no need to reset to zero)
      B.SetSubMatrix(0, 0, mgdshapex);
      B.SetSubMatrix(nd, 1, mgdshapey);
      B.SetSubMatrix(0, 2, mgdshapey);
      B.SetSubMatrix(nd, 2, mgdshapex);
#else
      // compute minimal tensor product needed for computation
      gxgxt = 0.;
      gygyt = 0.;
      gxgyt = 0.;
      gygxt = 0.;
      gxygxyt = 0.;
      AddMult_a_VVt(w, gdshapex, gxgxt);
      AddMult_a_VVt(w, gdshapey, gygyt);
      AddMult_a_VWt(w, gdshapex, gdshapey, gxgyt);
      gygxt.Transpose(gxgyt);
      gxygxyt = gxgyt;
      gxygxyt.Add(1., gygxt);
      /*
      cout << "gxgxt" << endl;
      gxgxt.Print(cout);
      cout << "gygyt" << endl;
      gygyt.Print(cout);
      cout << "gxgyt" << endl;
      gxgyt.Print(cout);
      cout << "gygxt" << endl;
      gygxt.Print(cout);
      cout << "gxygxyt" << endl;
      gxygxyt.Print(cout);
      */
#endif
      // get damage
      real_t d = dam.Eval(Tr, ip);
      // d=0.;

      // if non null damage
      if (d > 0.)
      {
         // to secure matrix non singular use limiter for damage
         d = min(d, 1. - limit);

         // computes gradiant of displacement
         MultAtB(disp, gdshape, strain);
         // cout << "disp" << endl;
         // disp.Print(cout);
         // cout << "disp^t.gdshape" << endl;
         // strain.Print(cout);
         // set as strain
         strain.Symmetrize();
         // cout << "strain" << endl;
         // strain.Print(cout);

#ifdef USE_AD
         vparam[0] = l;
         vparam[1] = m;
         vparam[2] = d;

         // cout<<ad_potentiel.Eval(vparam,vstrain)<<endl;
         ad_potentiel.Hessian(vparam, vstrain, hess);
         // cout<<"Hessian"<<endl;
         // hess.Print(cout);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) hook(i, j) = hess(i + 2 * (i % 2), j + 2 * (j % 2));
         hook(2, 2) = 0.5 * (hook(2, 2) + hess(1, 2));
         // cout<<"hook"<<endl;
         // hook.Print(cout);
#else

         // computing first and second invariant
         real_t I1 = strain(0, 0) + strain(1, 1);
         real_t I2 = strain(0, 1) * strain(0, 1) - strain(0, 0) * strain(1, 1);

         // if non null tensor
         if (I1 > limit || I2 > limit || I1 < mlimit || I2 < mlimit)
         {
            // computing the eigen values
            real_t delta = I1 * I1 + 4 * I2;
            MFEM_ASSERT(delta > -limit, "!!! Strain tensor not correctly defined: rank deficient !!!")
            real_t r = sqrt(max(0., delta));
            real_t e1 = (I1 + r) / 2.;
            real_t e2 = (I1 - r) / 2.;

            // cos and sin teta = delta_sigma/r and 2*sigma12/r, respectively
            real_t coss, sinn;
            if (r < limit)
            {  // singular second derivative for sigma_12->0 AND (sigma_11-sigma_22)->0
               real_t signe = (2 * strain(0, 1) / (strain(0, 0) - strain(1, 1))) > 0. ? 1 : -1;  // sign tangent
               coss = signe * sqrt(2.) / 2.;
               sinn = coss;
            }
            else
            {
               coss = (strain(0, 0) - strain(1, 1)) / r;
               sinn = 2 * strain(0, 1) / r;
            }
            // alpha's
            real_t alpha, alpha1, alpha2;
            if (e1 >= 0)
               alpha1 = 1.;
            else
               alpha1 = 0.;
            if (e2 >= 0)
               alpha2 = 1.;
            else
               alpha2 = 0.;
            if (I1 >= 0)
               alpha = 1.;
            else
               alpha = 0.;

            //  set constant coef
            real_t factor = 2. * m;      //  E / (1 + nu)
            real_t gamma = 0.5 * l / m;  // nu / (1 - 2 * nu)

            real_t c1 = 1. - alpha1 * d;
            real_t c2 = 1. - alpha2 * d;
            real_t c3 = 1. - alpha * d;

            // potentiel second derivative
            potentiel_sderivative(0, 0) = factor * (c1 + gamma * c3);
            potentiel_sderivative(0, 1) = factor * gamma * c3;
            potentiel_sderivative(1, 1) = factor * (c2 + gamma * c3);
            potentiel_sderivative(1, 0) = potentiel_sderivative(0, 1);
            // cout << "potentiel_sderivative" << endl;
            // potentiel_sderivative.Print(cout);

            // derivative eigenvalues
            dedeps(0, 0) = 1 + coss;
            dedeps(0, 1) = 1 - coss;
            dedeps(0, 2) = sinn;
            dedeps(1, 0) = 1 - coss;
            dedeps(1, 1) = 1 + coss;
            dedeps(1, 2) = -sinn;
            dedeps *= 0.5;
            // cout << "dedeps" << endl;
            // dedeps.Print(cout);

            // other terms
            real_t cos2 = coss * coss;
            real_t sin2 = sinn * sinn;
            real_t sc = sinn * coss;
            M(0, 0) = 1. - cos2;
            M(0, 1) = -1. + cos2;
            M(1, 1) = 1. - cos2;
            M(0, 2) = -1. * sc;
            M(1, 2) = 1. * sc;
            M(2, 2) = (1 - sin2);
            M(1, 0) = M(0, 1);
            M(2, 0) = M(0, 2);
            M(2, 1) = M(1, 2);
            M *= 0.5 * m;
            // cout << "M" << endl;
            // M.Print(cout);

            // tangent operator
            MultAtB(dedeps, potentiel_sderivative, hook);
            Mult(hook, dedeps, block);
            const real_t q = (r >= limit) ? (I1 / r * (c1 - c2) + (c1 + c2)) : (c1 + c2);
            Add(block, M, q, hook);
         }
         // null strain tensor
         else
         {
            // simple linear damaged hooke tensor
            hook = 0.;
            real_t md = (1. - d) * m;
            real_t ld = (1. - d) * l;
            hook(0, 0) = hook(1, 1) = 2 * md + ld;
            hook(1, 0) = hook(0, 1) = ld;
            hook(2, 2) = md;
         }
#endif
      }
      // if null damage
      else
      {
         // simple linear hook
         hook = 0.;
         hook(0, 0) = hook(1, 1) = 2 * m + l;
         hook(1, 0) = hook(0, 1) = l;
         hook(2, 2) = m;
      }
      // cout << "hook" << endl;
      // hook.Print(cout);

#ifdef USE_B
      Mult(B, hook, C);
      AddMult_a_ABt(w, C, B, elmat);
#else
      // check symmetry and null values
      MFEM_ASSERT(fabs(hook(2, 0) - hook(0, 2)) < limit, "non symmetric tensor ?? 2,0!=0,2")
      MFEM_ASSERT(fabs(hook(2, 1) - hook(1, 2)) < limit, "non symmetric tensor ?? 2,1!=1,2")
      MFEM_ASSERT(fabs(hook(1, 0) - hook(0, 1)) < limit, "non symmetric tensor ?? 1,0!=0,1")
      bool b20 = (fabs(hook(2, 0)) > limit);
      bool b21 = (fabs(hook(2, 1)) > limit);

      // xx diagonal block
      Add(hook(0, 0), gxgxt, hook(2, 2), gygyt, block);
      if (b20) Add(block, gxygxyt, hook(2, 0), block);
      elmat.AddSubMatrix(0, block);
      // yy diagonal block
      Add(hook(1, 1), gygyt, hook(2, 2), gxgxt, block);
      if (b21) Add(block, gxygxyt, hook(2, 1), block);
      elmat.AddSubMatrix(nd, nd, block);
      // xy off diagonal block
      Add(hook(1, 0), gxgyt, hook(2, 2), gygxt, block);
      if (b20) Add(block, gxgxt, hook(2, 0), block);
      if (b21) Add(block, gygyt, hook(1, 2), block);
      elmat.AddSubMatrix(0, nd, block);
      // yx off diagonal block
      block.Transpose();
      elmat.AddSubMatrix(nd, 0, block);
#endif
      // cout << "elmat" << endl;
      // elmat.Print(cout);
      time_grad += MPI_Wtime() - start;
   }
   void getTimes(double &v1, double &v2)
   {
      v1 += time_vect;
      v2 += time_grad;
   }

  private:
   Coefficient &lambda;
   Coefficient &mu;
   QuadratureFunctionCoefficient &dam;
   const IntegrationPoint &ip;
#ifndef USE_F_OUTSIDE
   VectorQuadratureFunctionCoefficient &load;
#endif
   double time_vect;
   double time_grad;
   const int nd;
   const int dim;
   real_t limit, mlimit;
#ifndef MFEM_THREAD_SAFE
   DenseMatrix dshape, gdshape, strain, sig, hook, block, disp, res;
#ifdef USE_B
   DenseMatrix B, C, mgdshapex, mgdshapey;
#else
   DenseMatrix gxgxt, gygyt, gxgyt, gygxt, gxygxyt;
#endif
#ifdef USE_AD
   DenseMatrix hess;
   Vector vparam;
   mfem::QFunctionAutoDiff<Potential, 4, 3> ad_potentiel;
#else
   DenseMatrix evec, potentiel_sderivative, dedeps, M;
   Vector eval, eigns;
#endif
   Vector vl, shape, f;
#endif
};

/// \brief test program
int main(int argc, char *argv[])
{
   // Initialize
   Mpi::Init(argc, argv);
   double measure[NB_MEASURE];
   std::fill(measure,measure+NB_MEASURE,0.);
   double start_all = MPI_Wtime();
   double start = MPI_Wtime();
   int nb_proc = Mpi::WorldSize();
   int mpi_rank = Mpi::WorldRank();

   // init output per proc
   string no = "proc_" + std::to_string(mpi_rank) + "_output.txt";
   if (mpi_rank > 50)
      freopen("/dev/null", "w", stdout);
   else
      freopen(no.c_str(), "w", stdout);

   const int order = 1;

   //== option treatement =========================================
#ifdef DEBUG_SQUARE
   const char *mesh_file = "data/square.msh";
#else
   const char *mesh_file = "data/neper_dam.msh";
#endif
   const char *petscrc_file = "";
   bool verbose = true;
   int max_refine = 0;

   OptionsParser args(argc, argv);
   // args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&petscrc_file, "-petscopts", "--petscopts", "PetscOptions file to use.");
   args.AddOption(&verbose, "-v", "--verbose", "-nv", "--not-verbose", "Output extra informations.");
   args.AddOption(&max_refine, "-r", "--refine", "Max refinement level.");
   real_t newton_rel_tol = 1e-7;
   real_t newton_abs_tol = 5e-8;
   args.AddOption(&newton_rel_tol, "-rel", "--relative-tolerance", "Relative tolerance for the Newton solve.");
   args.AddOption(&newton_abs_tol, "-abs", "--absolute-tolerance", "Absolute tolerance for the Newton solve.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi_rank == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi_rank == 0)
   {
      args.PrintOptions(cout);
      if (verbose) cout << "MFEM version " << GetVersionStr() << " (git sha:" << GetGitStr() << endl;
   }

#ifdef PETSC_SOLV
   MFEMInitializePetsc(NULL, NULL, petscrc_file, NULL);
#endif

   measure[1] += MPI_Wtime() - start;  // 1. Initialize
   start = MPI_Wtime();

   ParMesh *pmesh;
   int dim;
   {
      Mesh mesh(mesh_file, 1, 0, true);
      dim = mesh.Dimension();
      MFEM_ASSERT(dim == 2, "mesh file does not containe a 2D mesh")
      if (verbose)
      {
         cout << "Mesh nb tria:" << mesh.GetNE() << endl;
         cout << "Mesh nb edge:" << mesh.GetNEdges() << endl;
         cout << "Mesh nb node:" << mesh.GetNV() << endl;
         cout << "Mesh nb edge with ghost:" << mesh.GetNumFacesWithGhost() << endl;
         cout << "Mesh nb edge bnd:" << mesh.GetNFbyType(FaceType::Boundary) << endl;
         cout << "Mesh nb mat:" << mesh.attributes.Max() << endl;
         cout << "Mesh nb bnd:" << mesh.bdr_attributes.Max() << endl;
      }
      pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
   }
   measure[2] += MPI_Wtime() - start;  // 2.1 Read the mesh
   start = MPI_Wtime();
   // refine
   for (short i = 0; i < max_refine; ++i) pmesh->UniformRefinement(0);
   measure[3] += MPI_Wtime() - start;  // 2.2 Refine the mesh

   int nb_loc_elem = pmesh->GetNE();
   if (verbose)
   {
      cout << "Mesh nb tria loc:" << nb_loc_elem << endl;
      cout << "Mesh nb edge loc:" << pmesh->GetNEdges() << endl;
      cout << "Mesh nb node loc:" << pmesh->GetNV() << endl;
      cout << "Mesh nb edge with ghost local:" << pmesh->GetNumFacesWithGhost() << endl;
      cout << "Mesh nb edge bnd local:" << pmesh->GetNFbyType(FaceType::Boundary) << endl;
      cout << "Mesh nb mat loc:" << pmesh->attributes.Max() << endl;
      cout << "Mesh nb bnd loc:" << pmesh->bdr_attributes.Max() << endl;

#ifdef DEBUG_SQUARE
      AttributeSets &attr_sets = pmesh->attribute_sets;
      AttributeSets &bdr_attr_sets = pmesh->bdr_attribute_sets;
      std::set<string> names = attr_sets.GetAttributeSetNames();
      cout << "Element Attribute Set Names: ";
      for (auto const &set_name : names)
      {
         cout << " \"" << set_name << "\"";
      }
      cout << endl;
      std::set<string> bdr_names = bdr_attr_sets.GetAttributeSetNames();
      cout << "Boundary Attribute Set Names: ";
      for (auto const &bdr_set_name : bdr_names)
      {
         cout << " \"" << bdr_set_name << "\"";
      }
      cout << endl;
#endif
   }
   {
      start = MPI_Wtime();
      //== Constant poisson ==========================================
      real_t nu = 0.3;
      //== fill Young modulus ========================================
      std::vector<real_t> E_range(200);
      srand(6575);  // fix seed so that E_range is the same across process and for all runs
      for (auto &v : E_range)
      {
         const real_t a = (1.e8 - 5.e6) / 199.;
         v = a * (rand() % 200) + 5.e6;
#ifdef USE_ECST
         v = 1.e6;
#endif
      }
      //== fill lame as PWC ==========================================
      const real_t c1 = 1. + nu;
      const real_t c2 = nu / (c1 * (1. - 2. * nu));
      const real_t c3 = 1. / (2 * c1);
      int nb_mat = pmesh->attributes.Max();
      Vector lambda(nb_mat);
      Vector mu(nb_mat);
      for (int i = 0; i < nb_mat; ++i)
      {
         real_t Ec = E_range[(i + 1) % 200];
         lambda.Elem(i) = Ec * c2;
         mu.Elem(i) = Ec * c3;
      }
      PWConstCoefficient lambda_func(lambda);
      PWConstCoefficient mu_func(mu);
      measure[4] += MPI_Wtime() - start;  // 4.1 Material constant
      start = MPI_Wtime();
      //== Basic finite element ======================================
      //H1_FECollection fec(order, dim);
      LinearFECollection fec;
      //== Space =====================================================
      ParFiniteElementSpace space(pmesh, &fec, dim, Ordering::byVDIM);
      ParFiniteElementSpace spaceS(pmesh, &fec, 1);
      int nb_tv = spaceS.GetTrueVSize();
      int nb_v = spaceS.GetVSize();
      //== Quadrature Space ==========================================
      QuadratureSpace qspace1(pmesh, 1);
      QuadratureSpace qspace2(pmesh, 2);
      measure[5] += MPI_Wtime() - start;  // 3.1 Define space
      if (verbose)
      {
         auto nb_gtd = space.GlobalTrueVSize();
         auto nb_gd = space.GlobalVSize();
         if (mpi_rank == 0)
         {
            cout << "Number of true conforming dof: " << nb_gtd << endl;
            cout << "Number of dof: " << nb_gd << endl;
         }
         cout << "Number of local true conforming dof: " << space.GetTrueVSize() << endl;
         cout << "Number of local dof: " << space.GetVSize() << endl;
         cout << "Number of shared Faces:" << pmesh->GetNSharedFaces() << endl;
         // pmesh->bdr_attribute_sets.Print(cout);
         space.PrintPartitionStats();

         nb_gtd = spaceS.GlobalTrueVSize();
         nb_gd = spaceS.GlobalVSize();
         if (mpi_rank == 0)
         {
            cout << "Scalar : Number of true conforming dof: " << nb_gtd << endl;
            cout << "Scalar : Number of dof: " << nb_gd << endl;
         }
         cout << "Scalar : Number of local true conforming dof: " << nb_tv << endl;
         cout << "Scalar : Number of local dof: " << nb_v << endl;
         spaceS.PrintPartitionStats();
      }

      //== local dof id ==============================================
      ParGridFunction ids(&spaceS);
      Vector vv(nb_v);
      Vector vtv(nb_tv);
      for (int l = 0; l < nb_v; ++l) vv.Elem(l) = l;
      ids = vv;

      // mesh is made of same type of element thus the first one in local domain will give appropriate integration rule and
      // integration point for damIntegrator
      const IntegrationRule *ir1 = &qspace1.GetIntRule(0);
      const IntegrationPoint &ip1 = ir1->IntPoint(0);
      const IntegrationRule *ir2 = &qspace2.GetIntRule(0);

      start = MPI_Wtime();
      //== array to mark prop of interest ============================
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());

      //== Damage field construction =================================
      // identify dof with damge 1 using gmsh physical properties
      Array<int> dam_tdof_list;
      ess_bdr = 0;
#ifdef DEBUG_SQUARE
      // 4
      // 3,4,7
      ess_bdr[3] = 1;
#else
      // 148, 342, 333, 19, 380, 408, 328, 329, 325, 323, 96, 97, 531, 4,471, 234, 235, 184, 236, 419, 350, 332, 364,176,77,333,
      // 341, 343, 144, 143
      ess_bdr[147] = 1;
      ess_bdr[341] = 1;
      ess_bdr[332] = 1;
      ess_bdr[18] = 1;
      ess_bdr[379] = 1;
      ess_bdr[407] = 1;
      ess_bdr[327] = 1;
      ess_bdr[328] = 1;
      ess_bdr[324] = 1;
      ess_bdr[322] = 1;
      ess_bdr[95] = 1;
      ess_bdr[96] = 1;
      ess_bdr[530] = 1;
      ess_bdr[3] = 1;
      ess_bdr[470] = 1;
      ess_bdr[233] = 1;
      ess_bdr[234] = 1;
      ess_bdr[183] = 1;
      ess_bdr[235] = 1;
      ess_bdr[418] = 1;
      ess_bdr[349] = 1;
      ess_bdr[331] = 1;
      ess_bdr[363] = 1;
      ess_bdr[175] = 1;
      ess_bdr[76] = 1;
      ess_bdr[332] = 1;
      ess_bdr[340] = 1;
      ess_bdr[342] = 1;
      ess_bdr[143] = 1;
      ess_bdr[142] = 1;
#endif
      spaceS.GetEssentialTrueDofs(ess_bdr, dam_tdof_list);

      // create damage field using external data storage
      Vector v_dam_node(nb_v);
      ParGridFunction d(&spaceS, v_dam_node);

      // identify owned dof
      vtv = mpi_rank;
      d.SetFromTrueDofs(vtv);
      Vector v_rank(nb_v);
      v_rank = v_dam_node;

      // count edges
      Table &edge_vertex = *pmesh->GetEdgeVertexTable();
      Table vertex_edge;
      Transpose(edge_vertex, vertex_edge);
      const int nb_edges_loc = pmesh->GetNEdges();
      Array<int> shared(nb_edges_loc);
      shared = 0;
      const int nb_edges_shared = pmesh->GetNSharedFaces();
      for (int l = 0; l < nb_edges_shared; ++l) shared[pmesh->GetSharedFace(l)] = 1;
      vv = 0.;
      for (int l = 0; l < nb_v; ++l)
      {
         const int nb_edge = vertex_edge.RowSize(l);
         if (v_rank[l] == mpi_rank)
         {
            vv(l) = nb_edge;
         }
         else
         {
            const int *edges = vertex_edge.GetRow(l);
            int n_e = 0;
            for (int e = 0; e < nb_edge; ++e)
            {
               int e_id = edges[e];
               MFEM_ASSERT(edge_vertex.RowSize(e_id) == 2, "An edge got normaly 2 vertex !?")
               const int *vertex = edge_vertex.GetRow(e_id);
               int id = (vertex[0] == l) ? vertex[1] : vertex[0];
               if (v_rank[id] == mpi_rank || (shared[edges[e]] != 1)) ++n_e;
            }
            vv(l) = n_e;
         }
      }
      d = vv;
      Vector nb_edges_per_tnode(nb_tv);
      d.ParallelAssemble(nb_edges_per_tnode);
      // pre compute once the division
      for (int l = 0; l < nb_tv; ++l) nb_edges_per_tnode[l] = 1. / nb_edges_per_tnode[l];

      // set full damaged dof to MAX_DAM
      Vector v_dam(nb_tv);
      v_dam = 0.000;
      for (auto i : dam_tdof_list) v_dam(i) = MAX_DAM;
      d.SetFromTrueDofs(v_dam);

      // smooth
      for (short iter_smooth = 0, max_smooth = 8 * (max_refine + 1); iter_smooth < max_smooth; ++iter_smooth)
      {
         /*
         for (int l = 0; l < nb_tv; ++l) cout << l << " " << v_dam.Elem(l) << endl;
         cout << endl;
         for (int l = 0; l < nb_v; ++l) cout << l << " " << v_dam_node.Elem(l) << endl;
         cout << endl;
         */
         vv = 0.;
         for (int l = 0; l < nb_v; ++l)
         {
            real_t &di = v_dam_node.Elem(l);
            if (di < 0.01)
            {
               int nb_edge = vertex_edge.RowSize(l);
               // cout << l << " " << nb_edge << " " << di;
               real_t dam = 0.000;
               const int *edges = vertex_edge.GetRow(l);
               for (int e = 0; e < nb_edge; ++e)
               {
                  int e_id = edges[e];
                  const int *vertex = edge_vertex.GetRow(e_id);
                  int id = (vertex[0] == l) ? vertex[1] : vertex[0];
                  if (v_rank[id] == mpi_rank || (shared[edges[e]] != 1)) dam += v_dam_node.Elem(id);
                  // cout << " edge id " << e_id << " dofs " << vertex[0] << " " << vertex[1] << " dam " << dam;
               }
               vv[l] = dam;
               // cout << endl;
            }
         }
         d = vv;
         d.ParallelAssemble(vtv);
         for (int l = 0; l < nb_tv; ++l) v_dam[l] = std::max(vtv[l] * nb_edges_per_tnode[l], v_dam[l]);
         d.SetFromTrueDofs(v_dam);
         vv = 0.;
         for (int l = 0; l < nb_v; ++l)
         {
            int nb_edge = vertex_edge.RowSize(l);
            // real_t &di = v_dam_node.Elem(l);
            // cout << l << " " << nb_edge << " " << di;
            real_t dam = 0.000;
            const int *edges = vertex_edge.GetRow(l);
            for (int e = 0; e < nb_edge; ++e)
            {
               int e_id = edges[e];
               const int *vertex = edge_vertex.GetRow(e_id);
               int id = (vertex[0] == l) ? vertex[1] : vertex[0];
               if (v_rank[id] == mpi_rank || (shared[edges[e]] != 1)) dam += v_dam_node.Elem(id);
               // cout << " edge id " << e_id << " dofs " << vertex[0] << " " << vertex[1] << " dam " << dam;
            }
            vv[l] = dam;
            // cout << endl;
         }
         d = vv;
         d.ParallelAssemble(vtv);
         for (int l = 0; l < nb_tv; ++l) v_dam[l] = std::max(vtv[l] * nb_edges_per_tnode[l], v_dam[l]);
         d.SetFromTrueDofs(v_dam);
      }
      measure[6] += MPI_Wtime() - start;  // 4.2 Define damage
      // Compute at integration point d values once and encapsulate result into a coefficient via a quadrature function
      start = MPI_Wtime();
      QuadratureFunction qfd(qspace1, 1);
      qfd.ProjectGridFunction(d);
      QuadratureFunctionCoefficient qfcd(qfd);
      double qf = MPI_Wtime() - start;  // computing damage at ip is part of vector and matrix assembly
      measure[13] += qf / 2.;
      measure[14] += qf / 2.;
      start = MPI_Wtime();

      //== Solution on space =========================================
#if defined(IN_COMP) || defined(OUT_COMP)
      int nb_vxy = space.GetVSize();
      Vector v_x_node(nb_vxy);
      ParGridFunction x(&space, v_x_node);
#else
      ParGridFunction x(&space);
#endif

      //== use gmsh identified edges to impose Dirichlet =============
      // create all dirichlet dofs list
      Array<int> ess_tdof_list;
      ess_bdr = 0;
#ifdef DEBUG_SQUARE
      // left hand side
      ess_bdr[0] = 1;
      ess_bdr[4] = 1;
      // right hand side
      ess_bdr[1] = 1;
#else
      // left hand side
      ess_bdr[81] = 1;
      ess_bdr[188] = 1;
      ess_bdr[208] = 1;
      ess_bdr[263] = 1;
      ess_bdr[277] = 1;
      ess_bdr[484] = 1;
      ess_bdr[489] = 1;
      ess_bdr[503] = 1;
      ess_bdr[507] = 1;
      ess_bdr[533] = 1;
      ess_bdr[542] = 1;
      // right hand side
      ess_bdr[10] = 1;
      ess_bdr[49] = 1;
      ess_bdr[127] = 1;
      ess_bdr[204] = 1;
      ess_bdr[219] = 1;
      ess_bdr[457] = 1;
      ess_bdr[464] = 1;
      ess_bdr[488] = 1;
      ess_bdr[500] = 1;
      ess_bdr[524] = 1;
      ess_bdr[531] = 1;
      ess_bdr[536] = 1;
      ess_bdr[538] = 1;
#endif
      space.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

      // set non null Dirichlet associated to  x right end side component
      ess_bdr = 0;
#ifdef DEBUG_SQUARE
      ess_bdr[1] = 1;
#else
      ess_bdr[10] = 1;
      ess_bdr[49] = 1;
      ess_bdr[127] = 1;
      ess_bdr[204] = 1;
      ess_bdr[219] = 1;
      ess_bdr[457] = 1;
      ess_bdr[464] = 1;
      ess_bdr[488] = 1;
      ess_bdr[500] = 1;
      ess_bdr[524] = 1;
      ess_bdr[531] = 1;
      ess_bdr[536] = 1;
      ess_bdr[538] = 1;
#endif

#if 1
      // dirty  ??
      // get dof list
      Array<int> ess_tdof_list_2;
      // space.GetEssentialVDofs(ess_bdr, ess_tdof_list_2, 0);
      space.GetEssentialTrueDofs(ess_bdr, ess_tdof_list_2, 0);
      // create an intermediate vector
      Vector vd(space.GetTrueVSize());
      vd = 0.;
      // fill it for selected dof
      for (auto v : ess_tdof_list_2)
      {
#ifdef USE_TRAC
         vd[v] = 0.01;
#else
         vd[v] = -0.01;
#endif
      }
      // set x
      x.SetFromTrueDofs(vd);
#else
      // method by projection (more general)
#ifdef USE_TRAC
      ConstantCoefficient imp(0.01);
#else
      ConstantCoefficient imp(-0.01);
#endif
      Coefficient *coeff[1];
      coeff[0] = &imp;
      x.ProjectBdrCoefficient(coeff, ess_bdr);
#endif
      measure[7] += MPI_Wtime() - start;  // 5.2 Dirichlet setting
      start = MPI_Wtime();

      //== Arbitrary volumic load ====================================
      VectorFunctionCoefficient vl(dim, [](const Vector &x, Vector &f) {
         real_t xf = 100000.;
         // xf *= 1. / (0.0001 + x(0));
         real_t r = x(0) - 0.5;
         xf *= -r * r * r;
         // double xf = x(0);
         real_t y = x(1) - 0.5;
         f(0) = (1600. * y * y - 500.) * xf;
         f(1) = 0.;
      });
      // Use projection to be similar to FEniCSx interpolate
      // Using directly VectorFunctionCoefficient would be  more precise as true value would be used at
      // integration point but different from FEniCSx which use the linear interpolation from nodes value
      // But from cpu time point of view as load remain unchanged during NL iteration passing by
      // VectorQuadratureFunctionCoefficient would be more efficient a priory (computation at ip done only once)
      ParGridFunction load(&space);
      load.ProjectCoefficient(vl);
      measure[8] += MPI_Wtime() - start;  // 5.2 Neuman setting
      // computing load at ip have to be tracked
      start = MPI_Wtime();
      QuadratureFunction qfl(qspace2, 2);
      qfl.ProjectGridFunction(load);
      VectorQuadratureFunctionCoefficient qfcl(qfl);
      measure[13] += MPI_Wtime() - start;  // computing load at ip is part of vector assembly

#ifdef USE_F_OUTSIDE
#ifdef USE_VOLUME
      {
         start = MPI_Wtime();
         // use a linear form with conventionnel domain integrator to compute
         // volumic force
         ParLinearForm b(&space);
         auto integLoad = new VectorDomainLFIntegratorProf(qfcl);
         b.AddDomainIntegrator(integLoad);
         b.Assemble();
         // store in vd
         b.ParallelAssemble(vd);
         // reset to zero dirichlet so that they don't mess up the computation
         for (auto idx : ess_tdof_list) vd[idx] = 0.;
         // to be similar with FEniCSx count in 11 and 13 measure
         measure[11] += MPI_Wtime() - start;  // integrating loads is part of  nl resolution 
         measure[13] += integLoad->getMeasure();  // computing load vector is part of vector assembly count 
      }
#else
      vd.SetSize(0);
#endif
#else
      vd.SetSize(0);
#endif

      //== NonLinear form ============================================
      start = MPI_Wtime();
      ParNonlinearForm F(&space);
#ifdef USE_F_OUTSIDE
      auto pdi = new damIntegrator(lambda_func, mu_func, qfcd, ip1, ir2);
#else
      auto pdi = new damIntegrator(lambda_func, mu_func, qfcd, ip1, ir2, qfcl);
#endif
      F.AddDomainIntegrator(pdi);
      F.SetEssentialTrueDofs(ess_tdof_list);
      //== Vector on true dofs =======================================
      Vector X(space.GetTrueVSize());
      X = 0.;
      measure[9] += MPI_Wtime() - start;  // 7.1 Nonlinear form creation
      //== Linear solver ==============================================
      start = MPI_Wtime();
#ifdef PETSC_SOLV
      F.SetGradientType(Operator::Type::PETSC_MATAIJ);
      PetscPreconditioner prec(MPI_COMM_WORLD);
      PetscPCGSolver lin_solv(MPI_COMM_WORLD);
#else
      CGSolver lin_solv(MPI_COMM_WORLD);
      HypreBoomerAMG prec;
#ifdef PETSC_LIKE_BOOMERANG_SETTING
      // full manual setting via Hypre instance directely
      HYPRE_Solver prec_amg(prec);
      HYPRE_BoomerAMGSetNumFunctions(prec_amg, dim);
      HYPRE_BoomerAMGSetAggNumLevels(prec_amg, 0);        // no aggressive coarsening
      HYPRE_BoomerAMGSetStrongThreshold(prec_amg, 0.25);  //
      
      // full PETSC BOOMERANG setting (less efficiant and not keeped during all nl iter ?)
      //HYPRE_BoomerAMGSetCoarsenType(prec_amg, 6);         // Falgout
      //HYPRE_BoomerAMGSetInterpType(prec_amg, 0);          // classical
      //HYPRE_BoomerAMGSetRelaxType(prec_amg, 6);           // symmetric-SOR/Jacobi
      //HYPRE_BoomerAMGSetRelaxOrder(prec_amg, 1);          // Using CF-relaxation
                                                          
      prec.SetErrorMode(HypreSolver::ErrorMode::IGNORE_HYPRE_ERRORS);
#else
      prec.SetElasticityOptions(&space, true);
      // test
      // prec.SetStrengthThresh(0.25);
#endif
      // prec.SetPrintLevel(3);
#endif
      lin_solv.SetPreconditioner(prec);
      lin_solv.SetRelTol(1e-12);
      lin_solv.SetMaxIter(2000);
      lin_solv.SetPrintLevel(1);

      //== Non linear solver =========================================
      NewtonSolver nl_solv(MPI_COMM_WORLD);
      nl_solv.SetSolver(lin_solv);
      nl_solv.SetOperator(F);
      nl_solv.SetMaxIter(10);
      IterativeSolver::PrintLevel pl;
      pl.All();
      nl_solv.SetPrintLevel(pl);
      nl_solv.SetRelTol(newton_rel_tol);
      nl_solv.SetAbsTol(newton_abs_tol);
#ifdef TEST_ADAPT_LINRTOL
      nl_solv.SetAdaptiveLinRtol();
#endif
      measure[10] += MPI_Wtime() - start;  // 7.2 Solver creation
      start = MPI_Wtime();
      x.ParallelProject(X);
      nl_solv.Mult(vd, X);
      // MFEM_VERIFY(nl_solv.GetConverged(), "Newton Solver did not converge.");
      x.SetFromTrueDofs(X);
      measure[11] += MPI_Wtime() - start;  // 7.3 NonLinear resolution

      start = MPI_Wtime();
      //== strain computation ========================================
      DG_FECollection fecD3(0, 2);
      ParFiniteElementSpace spaceD3(
          pmesh, &fecD3, 3, Ordering::byVDIM);  // for 3 terms to store sym tensor: eps11 in x, eps12 in y, eps12=eps22 in z
      ParGridFunction strain(&spaceD3);
      strainTensor strain_coef(x);
      strain.ProjectCoefficient(strain_coef);
      //== stress computation ========================================
      ParGridFunction stress(&spaceD3);
      stressTensor stress_coef(x, qfcd, lambda_func, mu_func);
      stress.ProjectCoefficient(stress_coef);
      measure[15] += MPI_Wtime() - start;  // 8.1 strain/stress computation

      start = MPI_Wtime();
      //== Space for outputs =========================================
      DG_FECollection fecD(0, dim);
      ParFiniteElementSpace spaceD(pmesh, &fecD, 1);
      //== Partition for output ======================================
      ParGridFunction part(&spaceD);
      part = mpi_rank;

      //== Young modulus for output ==================================
      /* give only point with paraview
      QuadratureSpace qspace(pmesh,order);
      QuadratureFunction qfunc(qspace);
      Vector E(qspace.GetSize());
      for (int i = 0; i < nb_loc_elem; ++i)
      {
         E.Elem(i) = E_range[(pmesh->GetElement(i)->GetAttribute() - 1) % 200];
      }
      qfunc = E;
      */
      MFEM_ASSERT(nb_loc_elem == spaceD.GetVSize(), "Wrong discontinuous setting compared to nb_loc_elem")
      Vector E(spaceD.GetVSize());
      for (int i = 0; i < nb_loc_elem; ++i)
      {
         E.Elem(i) = E_range[(pmesh->GetElement(i)->GetAttribute()) % 200];
      }
      ParGridFunction Ef(&spaceD);
      Ef = E;

      //== Output ====================================================
      std::string out_filename;
      std::string sol_s = "nl_";
#ifdef USE_AD
      sol_s += "ad_";
#endif
#if defined(MFEM_USE_ADIOS2) && defined(USE_ADIOS_FOR_OUTPUT)
      // use adios
      // open stream
#ifdef USE_TRAC
      out_filename = "run/" + sol_s + "trac_" + std::to_string(nb_proc) + ".bp";
#else
      out_filename = "run/" + sol_s + "comp_" + std::to_string(nb_proc) + ".bp";
#endif
      adios2stream os(out_filename, adios2stream::openmode::out, MPI_COMM_WORLD, "BP4");
      //"BPFile"
      pmesh->Print(os);
      x.Save(os, "disp_mfem");
      part.Save(os, "part");
      Ef.Save(os, "E");
      d.Save(os, "d");
      load.Save(os, "load");
      strain.Save(os, "strain_mfem");
      stress.Save(os, "stress_mfem");
#else
      // use paraview
#ifdef USE_TRAC
      out_filename = "run/" + sol_s + "trac_" + std::to_string(nb_proc);
#else
      out_filename = "run/" + sol_s + "comp_" + std::to_string(nb_proc);
#endif
      {
         ParaViewDataCollection *pd = NULL;
         pd = new ParaViewDataCollection(out_filename, pmesh);
         pd->RegisterField("disp_mfem", &x);
         pd->RegisterField("part", &part);
         pd->RegisterField("E", &Ef);
         pd->RegisterField("d", &d);
         pd->RegisterField("ids", &ids);
         pd->RegisterField("load", &load);
         pd->RegisterField("strain_mfem", &strain);
         pd->RegisterField("stress_mfem", &stress);
         pd->SetLevelsOfDetail(order);
         pd->SetDataFormat(VTKFormat::BINARY);
         // pd->SetHighOrderOutput(true);
         pd->SetCycle(0);
         pd->SetTime(0.0);
         pd->Save();
         delete pd;
      }
#if 0
      // add native MFEM
#ifdef USE_TRAC
      out_filename = "run/"+sol_s+"trac_" + std::to_string(nb_proc) + ".sol";
#else
      out_filename = "run/"+sol_s+"comp_" + std::to_string(nb_proc) + ".sol";
#endif
      string mesh_ofile("run/mesh_" + std::to_string(nb_proc));
      ofstream mesh(mesh_ofile.c_str());
      pmesh->Print(mesh);
      ofstream sol(out_filename.c_str());
      x.Save(sol);
#endif
#endif
      measure[12] += MPI_Wtime() - start;     // 8 Outputs
      measure[0] += MPI_Wtime() - start_all;  // All

      pdi->getTimes(measure[13], measure[14]);

#ifdef OUT_COMP
      // for numerical comparison
      out_filename = "mfem_disp_mesh_" + std::to_string(nb_proc);
      ofstream ofile(out_filename.c_str());
      pmesh->Print(ofile);
      out_filename = "mfem_disp_sol_" + std::to_string(nb_proc);
      ofstream sol(out_filename.c_str());
      x.Save(sol);
#ifdef USE_VOLUME
      out_filename = "mfem_disp_" + std::to_string(max_refine) + "_vol";
#else
      out_filename = "mfem_disp_" + std::to_string(max_refine);
#endif
      ofstream oss(out_filename.c_str(), std::ofstream::out | std::ofstream::binary);
      int nbv = pmesh->GetNV();
      int q = 0;
      size_t siz = sizeof(real_t);
      for (int k = 0; k < nbv; ++k, q += 2)
      {
         real_t *coord = pmesh->GetVertex(k);
         oss.write(reinterpret_cast<const char *>(coord), siz);
         oss.write(reinterpret_cast<const char *>(coord + 1), siz);
         oss.write(reinterpret_cast<const char *>(&v_x_node[q]), siz);
         oss.write(reinterpret_cast<const char *>(&v_x_node[q + 1]), siz);
      }
      oss.close();
#endif
#ifdef IN_COMP
      //== field comparison with mfem non AD =========================
      {  // change scope to isolate this part
         assert(nb_proc == 1);
#ifdef USE_VOLUME
         string in_file_name = "data/mfem_disp_" + std::to_string(max_refine) + "_vol";
#else
         string in_file_name = "data/mfem_disp_" + std::to_string(max_refine);
#endif
         ifstream infile(in_file_name, std::ifstream::binary);
         real_t cx, cy, dx, dy;
         real_t errx = 0.;
         real_t erry = 0.;
         int q = 0;
         size_t siz = sizeof(real_t);
         // read raw data assuming that non AD version solution use the same mesh and that dof are rigorously in the same order
         infile.read(reinterpret_cast<char *>(&cx), siz);
         while (infile.good())
         {
            infile.read(reinterpret_cast<char *>(&cy), siz);
            infile.read(reinterpret_cast<char *>(&dx), siz);
            infile.read(reinterpret_cast<char *>(&dy), siz);
            real_t *coord = pmesh->GetVertex(q / 2);
            // cout<<q/2<<" "<<*coord<<" "<<cx<<" "<<*(coord+1)<<" "<<cy<<endl;
            assert(fabs(*coord - cx) < 1.e-6);
            assert(fabs(*(coord + 1) - cy) < 1.e-6);
            v_x_node[q] -= dx;
            v_x_node[q + 1] -= dy;
            // compute norme per component
            errx += v_x_node[q] * v_x_node[q];
            erry += v_x_node[q + 1] * v_x_node[q + 1];
            q += 2;
            infile.read(reinterpret_cast<char *>(&cx), siz);
         }
         assert(nb_vxy == q);
         real_t L2x = sqrt(errx);
         real_t L2y = sqrt(erry);
         cout << "Error L2 x:" << L2x << endl;
         cout << "Error L2 y:" << L2y << endl;
         ParGridFunction energy_error(&spaceD);
         energyError energy_coef(strain_coef, stress_coef);
         energy_error.ProjectCoefficient(energy_coef);
         cout << "Error in energy:" << energy_coef.getSum() << endl;
         energy_error *= 100. / energy_coef.getSum();
         L2x = 100. / L2x;
         L2y = 100. / L2y;
         q = 0;
         std::transform(v_x_node.begin(), v_x_node.end(), v_x_node.begin(), [&L2x, &L2y, &q](const double &v) {
            if ((q++) % 2) return v * L2y;
            return v * L2x;
         });
#ifdef USE_VOLUME
         out_filename = "run/diff_trac_vol" + std::to_string(max_refine);
#else
         out_filename = "run/diff_trac_vol" + std::to_string(max_refine);
#endif
         {
            ParaViewDataCollection *pd = NULL;
            pd = new ParaViewDataCollection(out_filename, pmesh);
            pd->RegisterField("AD-STD", &x);
            pd->RegisterField("Relative energy error AD-STD", &energy_error);
            pd->SetLevelsOfDetail(order);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetCycle(0);
            pd->SetTime(0.0);
            pd->Save();
            delete pd;
         }
      }
#endif
   }
   start_all = MPI_Wtime();
   delete pmesh;
   measure[0] += MPI_Wtime() - start_all;  // All
   double pow_measure[NB_MEASURE];
   double ecart_type_measure[NB_MEASURE];
   double max_measure[NB_MEASURE];
   double min_measure[NB_MEASURE];
   double avg_measure[NB_MEASURE];
   std::transform(measure, measure + NB_MEASURE, pow_measure, [](const double &v) { return v * v; });
   MPI_Reduce(measure, max_measure, NB_MEASURE, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(measure, min_measure, NB_MEASURE, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(measure, avg_measure, NB_MEASURE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pow_measure, ecart_type_measure, NB_MEASURE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (!mpi_rank)
   {
      cout << setfill('=') << setw(137) << " " << endl;
      cout << setfill(' ') << "| " << setw(12) << "val min"
           << " | " << setw(12) << "val max"
           << " | " << setw(12) << "std dev"
           << " | " << setw(12) << "% Coef var"
           << " | " << setw(12) << "val avg"
           << " | " << setw(12) << "% total avg"
           << " | " << setw(44) << " |" << endl;
      SL(0, "All")
      SL(1, "1. Initialize ")
      SL(2, "2.1 Read the mesh")
      SL(3, "2.2 Refine the mesh")
      SL(5, "3.1  Define space")
      SL(6, "3.2  Define damage")
      SL(4, "4.1  Material constant")
      SL(7, "5.1 Dirichlet setting")
      SL(8, "5.2 Neuman setting");
      SL(13, "6.1 Elementary vector");
      SL(14, "6.2 Elementary matrix");
      SL(9, "7.1 Nonlinear form creation");
      SL(10, "7.2 Solver creation");
      SL(11, "7.3 NonLinear resolution");
      SL(12, "8 Outputs");
      SL(15, "8.1 strain/stress computation");
      cout << setfill('=') << setw(137) << " " << endl;
   }

#ifdef PETSC_SOLV
   MFEMFinalizePetsc();
#endif
   Mpi::Finalize();
   return 0;
}
