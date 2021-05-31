surfaceScalarField Rhof
(
    "Rhof",
    linearInterpolate(rho)
);

surfaceVectorField Uf
(
    "Uf",
    linearInterpolate(U)
);

surfaceScalarField Pf
(
    "Pf",
    linearInterpolate(p)
);

surfaceScalarField Ef
(
    "Ef",
    linearInterpolate(e)
);

surfaceScalarField Hthermodynamicalf
(
    "Hthermodynamicalf",
    linearInterpolate(e)
);

surfaceScalarField MolWeightf
(
    "MolWeightf",
    linearInterpolate(MolWeight)
);

surfaceScalarField Inductionf
(
    "Inductionf",
    linearInterpolate(Induction)
);
//==========================================================================================
// lf - ninear face
surfaceScalarField Rholf
(
    "Rholf",
    linearInterpolate(rho)
);

surfaceVectorField Ulf
(
    "Ulf",
    linearInterpolate(U)
);

surfaceScalarField Plf
(
    "Plf",
    linearInterpolate(p)
);

surfaceScalarField Tlf
(
    "Tlf",
    linearInterpolate(T)
);

surfaceScalarField Elf
(
    "Elf",
    linearInterpolate(e)
);

surfaceScalarField MolWeightlf
(
    "MolWeightlf",
    linearInterpolate(MolWeight)
);

surfaceScalarField Inductionlf
(
    "Inductionlf",
    linearInterpolate(Induction)
);

surfaceScalarField Psilf
(
    "Psilf",
    linearInterpolate(psi)
);
//==========================================================================================
// Потоки уравнений неразрывности
//---------Start---------
surfaceScalarField RhoUf
(
    "RhoUf",
    (Rhof*Uf) & mesh.Sf()
);
//---------End---------


// Потоки уравнений импульса
//---------Start---------
surfaceVectorField RhoUUf
(
    "RhoUUf",
    (Rhof*Uf*Uf) & mesh.Sf()
);
surfaceVectorField PPf
(
    "PPf",
    Pf*mesh.Sf()
);
//---------End---------


// Потоки уравнений энергии
//---------Start---------
surfaceScalarField Hf
(
    "Hf",
    Ef + Pf/Rhof + 0.5*magSqr(Uf)
);
surfaceScalarField RhoUHf
(
    "HHf",
    (Rhof*Uf*Hf) & mesh.Sf()
);
//---------End---------


// Потоки уравнений молярной массы
//---------Start---------
surfaceScalarField RhoUMolWeightf
(
    "RhoUMolWeightf",
    (Rhof*Uf*MolWeightf) & mesh.Sf()
);
//---------End---------


// Потоки уравнений доли индукции
//---------Start---------
surfaceScalarField RhoUInductionf
(
    "RhoUInductionf",
    (Rhof*Uf*Inductionf) & mesh.Sf()
);
//---------End---------
