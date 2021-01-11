#ifndef HELMHOLTZCHANNEL_HPP
#define HELMHOLTZCHANNEL_HPP

//===============================
// HelmholtzSolver
//  - modifications:
//      1. initial condition to have turbulent velocity field
//      2. momentum source to drive flow (buoyancy-driven) 
//===============================

class HelmholtzBoussinesqSolver : public HelmholtzSolver { 
public:

  // MyHelmholtzSolver() {}
  HelmholtzBoussinesqSolver() {
  }

  void initData() { 
    HelmholtzSolver::initData();
  }

  // ~MyHelmholtzSolver() {}
  ~HelmholtzBoussinesqSolver() {
  }

  void initialHook() {
   
    if (step == 0) {

    if ( mpi_rank == 0 ) 
        cout << " >>>>> Calling initial hook" <<endl;

    // recall channel is -1 <= y <= 1...
    FOR_ICV {
      
    //  approximate turbulent mean profile...
 
     rho[icv] = 1.0;
//      rho[icv] = rho_const;
//    rho_uncorrelated[icv] = rho_const;
//    difference between rho and rho uncorrelated ?
    transport_scalar_vec[0][icv]=0.0;

    const double x  = x_cv[icv][0];
    const double absx = abs(x);

    const double ux = 10; 
    u[icv][0] = 0.0;
    u[icv][1] = 0.0;
    u[icv][2] = 0.0;

    u[icv][0] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);
    u[icv][1] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);
    u[icv][2] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);

// and rhoE...
//    rhoE[icv]  = p_ref/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
    }
  }

  }

  void temporalHook() {}

  void finalHook() {}

  HelmholtzBc* initHookBc(BfZone* p) {
    CERR(" > user must supply a hook boundary condition for zone: " << p->getName());
  }

  
  // the Helmholtz solver has implicit time advancement in a fractional 
  // step setting; as a result, the hooks for add source hooks are slightly
  // different.

  void momentumSourceHook(double * A,double (*rhs)[3]) {

    if ( mpi_rank == 0 ) 
      cout << ">>>>> adding momentum source, Boussinesq appriximation" << endl;

      const double T_ref = 0.0;
      const double beta = 0.0034; 
      const double g = 100;
      const double factor = 1.0;

//      transport_scalar_vec[0][icv]=50.0;
      FOR_ICV{
        rhs[icv][1] += factor*vol_cv[icv]*rho[icv]*g*beta*(transport_scalar_vec[0][icv]-T_ref);
      }

  }

  void massSourceHook(double * rhs) {}

};


#endif
