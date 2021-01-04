# CharLES Helmholtz solver + Boussinesq approximation (scalar tarnsport equation)
## 0. original solver + scalar boundary condition 

### validation case: turbulent channel flow
Q. Check governing equation?
    1. Schmidt number should go denominator (divided)?
    -> Varying trend of mean scalar profiles is the opposite to as it should be.

![scalar trend](channel_flow/images/scalar_profile_mean.png)

## 1. Add momentum source term (using Boussinesq approximation)

### To check code works: Rayleigh-Benard convection
![RB gif](Rayleigh-Benard/animation.gif)

### Validation case: Thermally-driven cavity


  
