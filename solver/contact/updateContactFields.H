// Описание функций
void Contact(Foam::rhoRDEThermo & thermo, const Foam::volScalarField & Induction, Foam::volScalarField & MolWeight, const Foam::volVectorField & U, Foam::volScalarField & Upotential,
    const volScalarField & Zmin, const Foam::volScalarField & Zmax,
    const volScalarField & Rmin, const Foam::volScalarField & Rmax,
    const Foam::volScalarField & InductionZmin, const Foam::volScalarField & InductionZmax,
    const Foam::volScalarField & InductionRmin, const Foam::volScalarField & InductionRmax,
    const Foam::volScalarField & FieldZmin, const Foam::volScalarField & FieldZmax,
    const Foam::volScalarField & FieldRmin, const Foam::volScalarField & FieldRmax,
    const Foam::volScalarField & FieldPhimin, const Foam::volScalarField & FieldPhimax,
    Foam::scalarField & ContactCoordinateFieldZ, const Foam::vectorField & ContactCoordinateMeshZ,
    Foam::scalarField & ContactCoordinateFieldR, const Foam::vectorField & ContactCoordinateMeshR,
    const Foam::scalar dt, const Foam::scalar CameraLength, const Foam::scalar CameraDiameter);

void UpdateContactField(const Foam::fvMesh & mesh,
    const Foam::volScalarField & Zmin, const Foam::volScalarField & Zmax, const Foam::volScalarField & Rmin, const Foam::volScalarField & Rmax, const Foam::volScalarField & Phimin, const Foam::volScalarField & Phimax,
    Foam::volScalarField & FieldZmin1, Foam::volScalarField & FieldZmax1, Foam::volScalarField & FieldZmin2, Foam::volScalarField & FieldZmax2,
    Foam::volScalarField & FieldRmin1, Foam::volScalarField & FieldRmax1, Foam::volScalarField & FieldRmin2, Foam::volScalarField & FieldRmax2,
    Foam::volScalarField & FieldPhimin, Foam::volScalarField & FieldPhimax,
    const Foam::surfaceScalarField & Field1f, const Foam::surfaceScalarField & Field2f, const Foam::surfaceScalarField & FieldPhif);
//==============================================================================
// Запуск функций


// Интерполяци по термодинамической части энергии
// Hthermodynamicalf -> Field
if(UseContact) UpdateContactField(mesh, Zmin, Zmax, Rmin, Rmax, Phimin, Phimax,    // Ячейки сетки
    InductionZmin, InductionZmax, FieldZmin, FieldZmax,     // Поле по оси Z
    InductionRmin, InductionRmax, FieldRmin, FieldRmax,     // Поле по оси R
    PPhimin, PPhimax,                                       // Поле по углу фи
    Inductionf, Hthermodynamicalf, Pf);  // Поля граней

// Интерполяция по молярной массе
// MolWeightf -> Field
//if(UseContact) UpdateContactField(mesh, Zmin, Zmax, Rmin, Rmax, Phimin, Phimax,    // Ячейки сетки
//    InductionZmin, InductionZmax, FieldZmin, FieldZmax,     // Поле по оси Z
//    InductionRmin, InductionZmax, FieldRmin, FieldRmax,     // Поле по оси R
//    PPhimin, PPhimax,                                       // Поле по углу фи
//    Inductionf, MolWeightf, Pf);  // Поля граней




if(UseContact && (runTime.value() > timeOnePerimeter/4.0) && !EmptyContactMesh)
    Contact(thermo, Induction, MolWeight, U, e,
        Zmin, Zmax, Rmin, Rmax, InductionZmin, InductionZmax, InductionRmin, InductionRmax,
        FieldZmin, FieldZmax, FieldRmin, FieldRmax, PPhimin, PPhimax,
        ContactCoordinateDataZ, ContactCoordinateMeshZ, ContactCoordinateDataR, ContactCoordinateMeshR,
        runTime.deltaTValue(), CameraLength, CameraDiameter);

