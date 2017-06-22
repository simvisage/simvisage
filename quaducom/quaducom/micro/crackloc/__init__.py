'''Package for the analysis of restrained localization.

module: mats1D_damage_x0.py

- One dimensional model of composite with strain softening behavior
  in a single material point at x = 0.
- The damage law has the standard shape
- The damage away from the x = 0 is given by the averaging function
  
  $\alpha = ( 1 - x^2 / R^2 )^2$  (1) 
  
  with R representing the interaction radius.
  Thus, the non-averaged damage is defined as 
  
  $\omega = \dirac(x) \omega_0$.    (2)
  
  The averaged value is then 
  
  \[
  \bar{\omega}(x) = \int_0^R \alpha(x,\xi) \cdot \omega d\xi
               = \alpha(x,0) \omega_0
  \]  (3)

- The package consists of one module specializing the mats1D_damage material 
  model. In particular, it introduces the
  * averaging function 
  * modifies the state update such that only the value for x = 0 is stored
  * other material points ($x \neq 0$) use the formula (3) to obtain their damage. 

- this model can be verified using a uniaxial bar with the damage returned as 
  trace.

module:
  
  crackloc.py

  define the idealization with the parameters for
  shape - fineness of discretization
  fe_type - fets for the idealization
  bond    - bond material model for the interaction between matrix and reinforcement.
  
'''