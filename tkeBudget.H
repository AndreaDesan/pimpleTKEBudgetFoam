dimensionedScalar currentTime("currentTime", dimTime, runTime.value());
dimensionedScalar currentTimeStep("currentTimeStep", dimTime, runTime.deltaTValue());

dimensionedScalar tAveraging("tAveraging", dimTime, currentTime.value()-tStartAveraging.value());

UAcc += U*currentTimeStep;
pAcc += p*currentTimeStep;

UTimeAveraged = UAcc/(max(tSmall,tAveraging));
pTimeAveraged = pAcc/(max(tSmall,tAveraging));

UPrime=(U-UTimeAveraged); //instantaneous velocity fluctuations vector
pPrime=(p-pTimeAveraged); //instantaneous pressure fluctuations

RSGS = turbulence->R(); //Instantaneous subgrid-scale reynolds stress tensor
kSGS = 0.5*tr(RSGS); //Instantaneous SGS turbulent kinetic energy
	
RRes = symm((UPrime*UPrime)); //Instantaneous resolved Reynolds stress tensor
kRes = 0.5*tr(RRes); //Instantaneous resolved tke

SRes = symm(fvc::grad(UPrime)); //Instantaneous resolved fluctuating rate of strain
epsRes=-2*nuLam*(SRes && SRes); //Instantaneous resolved dissipation rate

LESResIndex = kRes/max(kSmall,(kRes+kSGS)); //Instantaneous ratio of resolved to total (Resolved + SGS) tke

Info << "Statistics accumulated for " << tAveraging.value() << " s" << endl;
