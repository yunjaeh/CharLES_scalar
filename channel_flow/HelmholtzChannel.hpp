#ifndef HELMHOLTZCHANNEL_HPP
#define HELMHOLTZCHANNEL_HPP

//===============================
// HelmholtzSolver
//  - modifications:
//      1. initial condition to have turbulent velocity field
//      2. momentum source to compensate friction loss at walls
//===============================

class HelmholtzChannelSolver : public HelmholtzSolver { 
public:

  // MyHelmholtzSolver() {}
  HelmholtzChannelSolver() {}

  void initData() { 

    HelmholtzSolver::initData();
  }

  // ~MyHelmholtzSolver() {}
  ~HelmholtzChannelSolver() {}

  void initialHook() {
/*    
    if ( mpi_rank == 0 ) 
        cout << " >> Calling initial hook" <<endl;
   
    if (step == 0) {

    // recall channel is -1 <= y <= 1...
    FOR_ICV {
      
    //  approximate turbulent mean profile...
 
//   rho[icv] = 1.0;
      rho[icv] = rho_const;
//    rho_uncorrelated[icv] = rho_const;
//    difference between rho and rho uncorrelated ?

    const double y  = x_cv[icv][1];
    const double absy = abs(y);

//    const double ux = 18.0*(1.0 + y*y)*(1.0 - y*y)*(1.0 + y*y)*(1.0 - y*y); // ux = 18 in the center
    const double ux = -284.*pow(absy,7.)+356.*pow(absy,6.)-90.*pow(absy,4.)+20.; // ux = 20 in the center
  
    u[icv][0] = ux;
    u[icv][1] = 0.0;
    u[icv][2] = 0.0;

    u[icv][0] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);
    u[icv][1] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);
    u[icv][2] += ux*0.1*(double(rand())/double(RAND_MAX)-0.5);

// and rhoE...
//    rhoE[icv]  = p_ref/(gamma-1.0) + 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
    }
  }
*/
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
        cout << ">> adding momentum source" << endl;

    const double factor = 1.0;
    FOR_ICV{
        rhs[icv][0] += factor*vol_cv[icv]*rho[icv];
        }
    }
    void massSourceHook(double * rhs) {}

};


#endif
