# CharLES Helmholtz solver + Boussinesq approximation (scalar tarnsport equation)

### 0. original solver + wall boundary condition for scalar
##### validation case: ![turbulent channel flow](channel_flow/README.md)

Q. Check governing equation?
    1. Schmidt number should go denominator (divided)?
    -> Varying trend of mean scalar profiles is the opposite to as it should be.

![scalar trend](channel_flow/images/scalar_profile_mean.png)
    
    -> Add results using updated code


### 1. Add momentum source term (using Boussinesq approximation)
##### Test case to check code works: ![Rayleigh-Benard convection](Rayleigh-Benard/README.md)
![RB gif](Rayleigh-Benard/animation.gif)

##### Validation case: ![thermally-driven cavity](cavity/README.md)

### Comments
    1. Pressure and temperature coupling
    2. 
  
