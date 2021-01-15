CharLES Helmholtz solver + Boussinesq approximation 
(scalar transport equation)
----
- The objective is to check validity of using Helmholtz solver with scalar transport equation instead of ideal-gas solver, which solves for full temperature field (more expensive). 
- Boussinesq approximation is used to include momentum 
- First step is the validation of original solver (without buoyancy force), and the second step includes buoyancy effect

**Comments**

- Pressure and temperature coupling (their fluctuation)


## 1. Original solver + wall boundary condition for scalar

**validation case: [turbulent channel flow](channel_flow/)**

- Turbulent channel flow (*Re<sub>tau</sub>* =395)
- Constant temperature (scalar at upper and lower walls)
- No buoyancy effect

<img src="channel_flow/images/animation_U.gif" title="Development of streamwise velocity" width=400> <img src="channel_flow/images/animation_CT.gif" width=400>

<center>Figure: development of streamwise velocity & scalar field</center>



## 2. Add momentum source term (Boussinesq approximation)

**Test case to check code works: [Rayleigh-Benard convection](Rayleigh-Benard/)**

- Rayleigh-Benard convection: higher temperature at the bottom and lower temperature at the top boundary
- Instability occurs by temperature gradient in reverse to the gravity

<img src="Rayleigh-Benard/images/animation.gif" width=400>

<center>Figure: scalar field</center>



**Validation case: [thermally-driven cavity](cavity/)**

- Cavity flow driven by buoyancy force (temperature difference between walls)
- Run cases in different Rayleigh number (10<sup>4</sup>, 10<sup>5</sup>, 10<sup>6</sup>)

<img src="cavity/images/animation_U.gif" width=400><img src="cavity/images/animation_T.gif" width=400>

<center>Figure: velocity and temperature fields, Ra=10<sup>6</sup></center>



**Validation case: [Vertical turbulent channel flow](vertical_channel/)**

- Vertical channel frow driven by gravitational force (buoyancy) (temperature difference between walls)
- Run different Rayleigh numbers (Reference results available: DNS case)

<img src="vertical_channel/results/T_animation.gif" width=600>

<center>Figure: temperature field, Ra=5.4x10<sup>5</sup></center>



