// Функция отладочного отображения грани
void cellface(const Foam::fvMesh & mesh, Foam::label celli)
{
    
    const Foam::labelList & cFaces = mesh.cells()[celli];
    const Foam::point & centrevector = mesh.cellCentres()[celli];
    Foam::Info << "celli = " << celli << " - ";
    forAll(cFaces, cFacei)
    {
        Foam::label facei = cFaces[cFacei];
        const Foam::point & facevector = mesh.faceCentres()[facei];
        Foam::scalar delta = mag(facevector - centrevector);
        Info << facei << "|";
    }
    Foam::Info << endl;
}

template<typename TYPE>
inline void RenamePatchVolume(Foam::GeometricField<TYPE, Foam::fvPatchField, Foam::volMesh> & Field, word oldtype, word newtype)
{
    typename Foam::GeometricField<TYPE, Foam::fvPatchField, Foam::volMesh>::Boundary & Boundarys = Field.boundaryFieldRef();
    forAll(Boundarys, i)
    {
        if(Boundarys[i].type() == oldtype)
        {
            Foam::fvPatchField<TYPE> & PatchField = Boundarys[i];
            Boundarys.set(i, fvPatchField<TYPE>::New(newtype, oldtype, Field.mesh().boundary()[i], PatchField.internalField()));
        }
    }
};

template<typename TYPE>
inline void RenamePatchSurface(Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh> & Field, word oldtype, word newtype)
{
    typename Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Boundarys = Field.boundaryFieldRef();
    forAll(Boundarys, i)
    {
        if(Boundarys[i].type() == oldtype)
        {
            Foam::fvsPatchField<TYPE> & PatchField = Boundarys[i];
            Boundarys.set(i, fvsPatchField<TYPE>::New(newtype, oldtype, Field.mesh().boundary()[i], PatchField.internalField()));
        }
    }
};

// Искусственная жесткая стенка
void SolidWall(
const Foam::volVectorField & U, Foam::surfaceVectorField & Uf,
const Foam::volScalarField & P, Foam::surfaceScalarField & Pf,
const Foam::volScalarField & MolWeight, Foam::surfaceScalarField & MolWeightf,
const Foam::volScalarField & Induction, Foam::surfaceScalarField & Inductionf)
{
    const Foam::fvMesh & mesh = Uf.mesh();
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    const surfaceVectorField & FaceCentres = mesh.Cf();
    Foam::label deletecell = 0;
    forAll(owner, facei)
    {
        const Foam::point & FaceCentre = FaceCentres[facei];
        if(Foam::mag(FaceCentre.x()) < Epsilon && FaceCentre.y() > 0.0)
        {
            Foam::label cell = P[owner[facei]] < P[neighbour[facei]] ? owner[facei] : neighbour[facei];
            Pf[facei] = P[cell];
            MolWeightf[facei] = MolWeight[cell];
            Inductionf[facei] = Induction[cell];
            Uf[facei] = Foam::vector(0.0, 0.0, U[cell].z());
            deletecell = cell;
        }
    }

    forAll(mesh.boundaryMesh(), ipatch)
    {
        Foam::word BCtype = mesh.boundaryMesh().types()[ipatch];
        Foam::fvsPatchVectorField & Usf = Uf.boundaryFieldRef()[ipatch];
        Foam::fvsPatchScalarField & Psf = Pf.boundaryFieldRef()[ipatch];
        Foam::fvsPatchScalarField & MolWeightsf = MolWeightf.boundaryFieldRef()[ipatch];
        Foam::fvsPatchScalarField & Inductionsf = Inductionf.boundaryFieldRef()[ipatch];

        if (mesh.boundaryMesh()[ipatch].coupled())
        {
            //Foam::tmp<Foam::scalarField> PInternal = P.boundaryField()[ipatch].patchInternalField();
            //Foam::tmp<Foam::scalarField> PNeighbour = P.boundaryField()[ipatch].patchNeighbourField();
            //Foam::tmp<Foam::scalarField> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            //Foam::tmp<Foam::scalarField> MolWeightNeighbour = MolWeight.boundaryField()[ipatch].patchNeighbourField();
            //Foam::tmp<Foam::scalarField> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            //Foam::tmp<Foam::scalarField> InductionNeighbour = Induction.boundaryField()[ipatch].patchNeighbourField();
            const Foam::fvsPatchVectorField & FaceCentres = mesh.Cf().boundaryField()[ipatch];
            forAll(FaceCentres, facei)
            {
                const Foam::point & FaceCentre = FaceCentres[facei];
                if(Foam::mag(FaceCentre.x()) < Epsilon*100.0 && FaceCentre.y() > 0.0)
                {
                    //Psf[facei] = MIN(PInternal.ref()[facei], PNeighbour.ref()[facei]);
                    //MolWeightsf[facei] = MIN(MolWeightInternal.ref()[facei], MolWeightNeighbour.ref()[facei]);
                    //Inductionsf[facei] = MAX(InductionInternal.ref()[facei], InductionNeighbour.ref()[facei]);
                    Usf[facei] = Foam::vector(0.0, 0.0, Usf[facei].z());
                    Psf[facei] = P[deletecell];
                    MolWeightsf[facei] = MolWeight[deletecell];
                    Inductionsf[facei] = Induction[deletecell];
                }
            }
        }else
        {
            const Foam::fvsPatchVectorField & FaceCentres = mesh.Cf().boundaryField()[ipatch];
            forAll(FaceCentres, facei)
            {
                const Foam::point & FaceCentre = FaceCentres[facei];
                if(Foam::mag(FaceCentre.x()) < Epsilon*100 && FaceCentre.y() > 0.0)
                {
                    Foam::label cell = facei;
                    Usf[facei] = Foam::vector(0.0, 0.0, U.boundaryField()[ipatch][cell].z());
                    Psf[facei] = P.boundaryField()[ipatch][cell];
                    MolWeightsf[facei] = MolWeight.boundaryField()[ipatch][cell];
                    Inductionsf[facei] = Induction.boundaryField()[ipatch][cell];
                }
            }
        }
        
    }
}

//=========================================================================================================
// Последний вход в искусственную жесткую стенку
void LastSolidWall(Foam::volScalarField & E, Foam::volScalarField & T, Foam::volVectorField & U, Foam::volScalarField & P, Foam::volScalarField & MolWeight, Foam::volScalarField & Induction)
{
    Foam::Info << "Disappear Solid Wall" << Foam::endl;
    Foam::Info << "***" << Foam::endl;
    const Foam::fvMesh & mesh = U.mesh();
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    const surfaceVectorField & FaceCentres = mesh.Cf();
    Foam::label Plarge = 0;
    Foam::label Psmall = 0;
    forAll(owner, facei)
    {
        const Foam::point & FaceCentre = FaceCentres[facei];
        if(Foam::mag(FaceCentre.x()) < Epsilon && FaceCentre.y() > 0.0)
        {
            Plarge = P[owner[facei]] > P[neighbour[facei]] ? owner[facei] : neighbour[facei];
            Psmall = P[owner[facei]] < P[neighbour[facei]] ? owner[facei] : neighbour[facei];
            // Ячейки
            U[Plarge] = U[Psmall];
            P[Plarge] = P[Psmall];
            MolWeight[Plarge] = MolWeight[Psmall];
            Induction[Plarge] = Induction[Psmall];
            T[Plarge] = T[Psmall];
            E[Plarge] = E[Psmall];
        }
    }

    // Граница
    forAll(mesh.boundaryMesh(), ipatch)
    {
        Foam::word BCtype = mesh.boundaryMesh().types()[ipatch];
        if (mesh.boundaryMesh()[ipatch].coupled())
        {
            const Foam::fvsPatchVectorField & FaceCentres = mesh.Cf().boundaryField()[ipatch];
            //----------------------------------------------------------
            Foam::tmp<Foam::scalarField> tmpEInternal(E.boundaryFieldRef()[ipatch].patchInternalField());
            Foam::tmp<Foam::scalarField> tmpPInternal(P.boundaryFieldRef()[ipatch].patchInternalField());
            Foam::tmp<Foam::scalarField> tmpTInternal(T.boundaryFieldRef()[ipatch].patchInternalField());
            Foam::tmp<Foam::vectorField> tmpUInternal(U.boundaryFieldRef()[ipatch].patchInternalField());
            Foam::tmp<Foam::scalarField> tmpMolWeightInternal(MolWeight.boundaryFieldRef()[ipatch].patchInternalField());
            Foam::tmp<Foam::scalarField> tmpInductionInternal(Induction.boundaryFieldRef()[ipatch].patchInternalField());
            //----------------------------------------------------------
            Foam::scalarField & EInternal = tmpEInternal.ref();
            Foam::scalarField & PInternal = tmpPInternal.ref();
            Foam::scalarField & TInternal = tmpTInternal.ref();
            Foam::vectorField & UInternal = tmpUInternal.ref();
            Foam::scalarField & MolWeightInternal = tmpMolWeightInternal.ref();
            Foam::scalarField & InductionInternal = tmpInductionInternal.ref();
            forAll(FaceCentres, facei)
            {
                const Foam::point & FaceCentre = FaceCentres[facei];
                if(Foam::mag(FaceCentre.x()) < Epsilon && FaceCentre.y() > 0.0)
                {
                    EInternal[facei] = E[Psmall];
                    PInternal[facei] = P[Psmall];
                    TInternal[facei] = T[Psmall];
                    UInternal[facei] = U[Psmall];
                    MolWeightInternal[facei] = MolWeight[Psmall];
                    InductionInternal[facei] = Induction[Psmall];
                }
            }
            E.boundaryFieldRef()[ipatch].patchInternalField(EInternal);
            P.boundaryFieldRef()[ipatch].patchInternalField(PInternal);
            T.boundaryFieldRef()[ipatch].patchInternalField(TInternal);
            U.boundaryFieldRef()[ipatch].patchInternalField(UInternal);
            MolWeight.boundaryFieldRef()[ipatch].patchInternalField(MolWeightInternal);
            Induction.boundaryFieldRef()[ipatch].patchInternalField(InductionInternal);
        }
    }
}
