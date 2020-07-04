# pimpleTKEBudgetFoam
Extension of the OpenFOAM's pimpleFoam solver to calculate the total (SGS + resolved) turbulent kinetic energy ![k](http://latex.codecogs.com/gif.latex?k) and turbulent dissipation rate ![\varepsilon](https://latex.codecogs.com/gif.latex?\varepsilon) at runtime during LES simulations. It can be easily extended to include all the terms of the turbulent kinetic energy budget.

In addition, the LES resolution index

![LES_{ResIndex}=\frac{k_{Res}}{k_{Res}+k_{SGS}}](http://latex.codecogs.com/gif.latex?LES_%7BResIndex%7D%3D%5Cfrac%7Bk_%7BRes%7D%7D%7Bk_%7BRes%7D&plus;k_%7BSGS%7D%7D)

is also calculated at runtime.

Please note that the utlity assumes that the fieldAverage utility is used to calculate UMean. This is then used to calculate the fluctuating velocity vector UPrime within tkeBudget.H as UPrime = U - UMean.

If the time-averaging of U via fieldAverage is not set up in controlDict, the solver raises a warning and both kTot and epsilonTot will not be updated at runtime, keeping their initial value of zero (see the intialization in createFields.H and the "if" statement in tkeBudget.H).
