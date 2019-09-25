# pimpleTKEBudgetFoam
Extension of the OpenFOAM's pimpleFoam solver to calculate the resolved Reynolds stress tensor, turbulent kinetic energy ![k](http://latex.codecogs.com/gif.latex?k) and turbulent dissipation rate ![$\varepsilon$](http://latex.codecogs.com/gif.latex?%24%5Cvarepsilon%24) at runtime during LES simulations. It can be easily extended to include all the terms of the turbulent kinetic energy budget.

In addition, the LES resolution index

![LES_{ResIndex}=\frac{k_{Res}}{k_{Res}+k_{SGS}}](http://latex.codecogs.com/gif.latex?LES_%7BResIndex%7D%3D%5Cfrac%7Bk_%7BRes%7D%7D%7Bk_%7BRes%7D&plus;k_%7BSGS%7D%7D)

is also calculated at runtime.

Please note that:
 - The time-averaging of the instantaneous velocity and pressure is hardcoded in the solver. The resulting time-averaged fields are UTimeAveraged and pTimeAveraged
 - The averaging of other desired quantities has to be set up via the fieldAverage utility
 - To allow for the averaging to span over different runs, the starting time of the averaging process has to be specified via the tStartAveraging entry in the averagingProperties dictionary (an example of the dictionary is included in the repository)
