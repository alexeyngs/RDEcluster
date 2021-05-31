Info << "Reading thermophysical properties\n" << endl;

autoPtr<Foam::rhoRDEThermo> pThermo
(
    Foam::rhoRDEThermo::New(mesh)
);

Foam::rhoRDEThermo & thermo = pThermo();

thermo.correct();

Foam::volScalarField & e = thermo.he();
Foam::volScalarField & p = thermo.p();
Foam::volScalarField & T = thermo.T();
const Foam::volScalarField & psi = thermo.psi();
const Foam::volScalarField & mu = thermo.mu();
Foam::volScalarField & MolWeight = thermo.MolWeight();
Foam::volScalarField & FMolWeight = thermo.FMolWeight();
Foam::volScalarField & Induction = thermo.Induction();
Foam::volScalarField & FInduction = thermo.FInduction();

Foam::volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Foam::volVectorField U
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

Foam::volVectorField Ucylinder
(
    IOobject
    (
        "Ucylinder",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);

// OpenFOAM delta coeffs can be negative for coupled patches
Foam::volScalarField hdelta
(
    IOobject
    (
        "hdelta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1.0
);

Foam::scalar MinHdelta = Foam::GREAT;
forAll(mesh.cells(), celli) // По всем ячейкам
{
    Foam::scalar lenght = Foam::GREAT; // Максимально возможный размер
    const Foam::labelList & cFaces = mesh.cells()[celli];
    const Foam::point & centrevector = mesh.cellCentres()[celli];
    forAll(cFaces, cFacei)
    {
        Foam::label facei = cFaces[cFacei];
        const Foam::point & facevector = mesh.faceCentres()[facei];
        Foam::scalar delta = mag(facevector - centrevector);
        if(lenght > delta) lenght = delta;
    }
    hdelta[celli] = lenght;
    if(MinHdelta > lenght) MinHdelta = lenght;
}

//------------------------------
void MultiplyWave(const Foam::label Multiplier, const Foam::scalar CameraLength, const Foam::scalar MinHdelta, const Foam::fvMesh & mesh, Foam::rhoRDEThermo & thermo, Foam::volScalarField & rho, Foam::volVectorField & U);
if(Multiplier > 1) MultiplyWave(Multiplier, CameraLength, MinHdelta, mesh, thermo, rho, U);
//------------------------------

Foam::volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

Foam::volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

Foam::volScalarField rhoMolWeight
(
    IOobject
    (
        "rhoMolWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*MolWeight
);

Foam::volScalarField rhoInduction
(
    IOobject
    (
        "rhoInduction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Induction
);


Foam::volScalarField gamma
(
    "gamma",
    thermo.Cp()/thermo.Cv()
);

Foam::volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(gamma  / psi)
);

Foam::volScalarField Mz
(
    IOobject
    (
        "Mz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Z) / c
);

Foam::volScalarField Courant
(
    "Courant",
    (mag(U) + c)*runTime.deltaT() / hdelta
);
//==========================================================================================================

Foam::volScalarField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U.component(vector::Z)
);

Foam::volScalarField impuls
(
    IOobject
    (
        "impuls",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p + rho*U.component(vector::Z)*U.component(vector::Z)
);
