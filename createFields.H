Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const dimensionedScalar tStart("tStart", runTime.startTime());
const dimensionedScalar tSmall("tSmall", dimTime, small);

//user must specify the time at which sampling started in the averagingPropertis dictionary
IOdictionary averagingProperties 
(
    IOobject
    (
        "averagingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar tStartAveraging
(
    "tStartAveraging",
    dimTime,
    averagingProperties.lookup("tStartAveraging")
);

Info << "Started to collect statistics at t = " << tStartAveraging.value() << " s" << endl;

//TimeAveraged fields are read if present 

volVectorField UTimeAveraged
(
    IOobject
    (
        "UTimeAveraged",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "UTimeAveraged", 
           dimVelocity, 
           vector (0,0,0)
      )
);

volScalarField pTimeAveraged
(
    IOobject
    (
        "pTimeAveraged",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedScalar
      (
           "pTimeAveraged", 
           dimVelocity*dimVelocity, 
           scalar(0.0)
      )    

);

//Instantiate accumulated U*time_step and p*time_step to calculate the average
//Initialize as Averaged_Value/Sampling_time (where sampling time = tStart - tStartAveraging) to continue averaging from previous fields
volVectorField UAcc
(
    IOobject
    (
        "UAcc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        UTimeAveraged*(tStart-tStartAveraging)
);

volScalarField pAcc
(
    IOobject
    (
        "pAcc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        pTimeAveraged*(tStart-tStartAveraging)
);

//Instantiate instantaneous values of velocity and pressure fluctuations
volVectorField UPrime
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        U-UTimeAveraged
);

volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        p-pTimeAveraged
);



#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


//Instantiate useful instantaneous fields to calculate the budgets
//The corresponding time-averaged quantities for the various fields can be calculated using the fieldAverage function object
const volScalarField nuLam(turbulence->nu());

const volScalarField kSmall
(
    IOobject
    (
        "kSmall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
           dimensionedScalar
      (
           "kSmall", 
           dimVelocity*dimVelocity,
           SMALL 
      )
);

volSymmTensorField RSGS
(
    IOobject
    (
        "RSGS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        turbulence->R()
);


volScalarField kSGS
(
    IOobject
    (
        "kSGS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        0.5*tr(RSGS)
);

volSymmTensorField RRes
(
    IOobject
    (
        "RRes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        symm((UPrime*UPrime))
);

volScalarField kRes
(
    IOobject
    (
        "kRes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5*tr(RRes)
);

volSymmTensorField SRes
(
    IOobject
    (
        "Sres",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        symm(fvc::grad(UPrime))
);

volScalarField epsRes
(
    IOobject
    (
        "epsRes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        -2*nuLam*(SRes && SRes)
);

volScalarField LESResIndex 
(
    IOobject
    (
        "LESResIndex",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        kRes/max(kSmall,(kRes+kSGS))
);

#include "createMRF.H"
#include "createFvOptions.H"
