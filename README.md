# CharLES Helmholtz solver + Boussinesq approximation (scalar tarnsport equation)

## 0. original solver + wall boundary condition for scalar
### validation case: [turbulent channel flow](channel_flow/)

1. Schmidt number should go denominator (divided)? \
-> Varying trend of mean scalar profiles is the opposite to as it should be.
<img src="channel_flow/images/scalar_profile_mean.png" width=600>

Updated results:


code modification:
    **HelmholtzSolver.cpp & HelmholtzSolverBCs.cpp **
    Modify scalar transport equation part:
    Schmidt number is divided instead of being multiplied     

    **FlowSolver.hpp**
    Initial value of Schmidt number from 0.0 to 1.0 \
    to avoid errors when a scalar is initialized


## 1. Add momentum source term (using Boussinesq approximation)
### Test case to check code works: [Rayleigh-Benard convection](Rayleigh-Benard/README.md)
1) Rayleigh-Benard convection: higher temperature at the bottom and lower temperature at the top boundary
2) Instability occurs by temperature gradient

<img src="Rayleigh-Benard/animation.gif" width=600>


##### Validation case: [thermally-driven cavity](cavity/)
1) Cavity flow driven by buoyancy force (temperature difference between walls)
2) Run cases in different Rayleigh number (10^4, 10^5, 10^6)

<img src="cavity/images/results_Ra_10_6.png" width=600>


### Comments
1. Pressure and temperature coupling (their fluctuation)
2. 
  
