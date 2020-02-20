if (U.db().foundObject<volVectorField>("UMean"))
{

    const objectRegistry& db = U.db();

    volVectorField UPrime = (U - db.lookupObject<volVectorField>("UMean"));

    volScalarField kSGS = turbulence->k(); //Instantaneous subgrid-scale tke
	    
    volScalarField kRes = 0.5*(UPrime & UPrime); //Instantaneous resolved tke

    kTot = kRes + kSGS;

    volScalarField epsilonSGS = turbulence->epsilon();

    volSymmTensorField SRes = symm(fvc::grad(UPrime)); //Instantaneous resolved fluctuating rate of strain
    volScalarField epsilonRes=-2*nuLam*(SRes && SRes); //Instantaneous resolved dissipation rate

    epsilonTot = mag(epsilonRes + epsilonSGS);

    LESResIndex = kRes/max(kSmall,kTot); //Instantaneous ratio of resolved to total (Resolved + SGS) tke



    Info << "ESTIMATED LESResIndex min, max, avg: " << min(LESResIndex).value()
					    << ", " << max(LESResIndex).value()
					    << ", " << LESResIndex.weightedAverage(mesh.V()).value()
					    << endl;

}

else
{
     Warning << "Averaging of velocity U is not active in controlDict." << endl 
             << "This is necessary to calculate resolved TKE and dissipation rate." << endl
             << "Add field U to the fieldAverage utility iin controlDict." << endl;
}
