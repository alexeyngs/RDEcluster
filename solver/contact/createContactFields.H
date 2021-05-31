
Foam::volScalarField Zmin
(
    IOobject
    (
        "Zmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    CameraLength
);

Foam::volScalarField Zmax
(
    IOobject
    (
        "Zmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE 
    ),
    mesh,
    0.0
);

Foam::volScalarField Rmin
(
    IOobject
    (
        "Rmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    CameraDiameter * 0.5
);

Foam::volScalarField Rmax
(
    IOobject
    (
        "Rmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE 
    ),
    mesh,
    0.0
);

Foam::volScalarField Phimin
(
    IOobject
    (
        "Phimin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    2.0*Foam::constant::mathematical::pi
);

Foam::volScalarField Phimax
(
    IOobject
    (
        "Phimax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);
//-------------------------------------------------------------------------------------

Foam::volScalarField InductionZmin
(
    IOobject
    (
        "InductionZmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField InductionZmax
(
    IOobject
    (
        "InductionZmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField InductionRmin
(
    IOobject
    (
        "InductionRmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField InductionRmax
(
    IOobject
    (
        "InductionRmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldZmin
(
    IOobject
    (
        "FieldZmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldZmax
(
    IOobject
    (
        "FieldZmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldRmin
(
    IOobject
    (
        "FieldRmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldRmax
(
    IOobject
    (
        "FieldRmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField PPhimin
(
    IOobject
    (
        "PPhimin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField PPhimax
(
    IOobject
    (
        "PPhimax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);
// Описание функций
void GetContactEdgeSize(const Foam::fvMesh & mesh, const Foam::scalar CameraLength, const Foam::scalar CameraDiameter, bool & EmptyContactMesh,
    Foam::vectorField & ContactCoordinateMeshZ, Foam::scalarField & ContactCoordinateFieldZ,
    Foam::vectorField & ContactCoordinateMeshR, Foam::scalarField & ContactCoordinateFieldR);



void UpdateContactMesh(
    Foam::volScalarField & Zmin, Foam::volScalarField & Zmax,
    Foam::volScalarField & Rmin, Foam::volScalarField & Rmax,
    Foam::volScalarField & Phimin, Foam::volScalarField & Phimax);


//---------------------------------------------------------------------------------------------------------------------------
// Поля

// Z
Foam::vectorField ContactCoordinateMeshZ;    // Сетка с координатами
Foam::scalarField ContactCoordinateDataZ;    // Сетка с данными
// R
Foam::vectorField ContactCoordinateMeshR;    // Сетка с координатами
Foam::scalarField ContactCoordinateDataR;    // Сетка с данными
// Вызов функций - получение координат и данных
GetContactEdgeSize(mesh, CameraLength, CameraDiameter, EmptyContactMesh,
    ContactCoordinateMeshZ, ContactCoordinateDataZ, ContactCoordinateMeshR, ContactCoordinateDataR);

UpdateContactMesh(Zmin, Zmax, Rmin, Rmax, Phimin, Phimax);

