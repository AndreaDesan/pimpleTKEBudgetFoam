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
volScalarField kTot
(
    IOobject
    (
        "kTot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
     mesh,
           dimensionedScalar
      (
           "kTot",
           dimVelocity*dimVelocity,
           0.0
      )
);

volScalarField epsilonTot
(
    IOobject
    (
        "epsilonTot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
     mesh,
           dimensionedScalar
      (
           "epsilonTot",
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

volScalarField LESResIndex
(
    IOobject
    (
        "LESResIndex",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
     mesh,
           dimensionedScalar
      (
           "LESResIndex",
           dimless,
           0.0
      )
);

#include "createMRF.H"
#include "createFvOptions.H"
