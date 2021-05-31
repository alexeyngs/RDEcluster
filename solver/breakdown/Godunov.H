template<typename TScalar, typename TVector, typename SFType>
inline void Godunov
(
    Foam::rhoRDEThermo & Thermo,
    const SFType & Sf,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & Rho,
    const Foam::GeometricField<TVector, Foam::fvPatchField, Foam::volMesh> & U,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & p,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & E,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & MolWeight,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & Induction,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & c,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & gamma,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & T,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rhof,
    Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Uf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Pf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ef,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Hthermodynamicalf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionf
)
{
    const Foam::fvMesh & mesh = Rho.mesh();
    const Foam::labelUList & P = mesh.owner();
    const Foam::labelUList & N = mesh.neighbour();
    //-------------------------------------------------------------------------------
    const Foam::tmp<surfaceScalarField> & tlambdas = mesh.surfaceInterpolation::weights();
    if (surfaceInterpolation::debug)
    {
        InfoInFunction << "Interpolating " << Rho.type() << " " << Rho.name() << " from cells to faces without explicit correction" << endl;
    }

    const Foam::surfaceScalarField & lambdas = tlambdas();
    const Foam::scalarField & lambda = lambdas;
    
    const Foam::Field<TScalar> & Rho_ = Rho;
    const Foam::Field<TVector> & U_ = U;
    const Foam::Field<TScalar> & p_ = p;
    const Foam::Field<TScalar> & E_ = E;
    const Foam::Field<TScalar> & MolWeight_ = MolWeight;
    const Foam::Field<TScalar> & Induction_ = Induction;
    const Foam::Field<TScalar> & c_ = c;
    const Foam::Field<TScalar> & gamma_ = gamma;
    const Foam::Field<TScalar> & T_ = T;
    
    
    Foam::Field<TScalar> & Rhof_ = Rhof.primitiveFieldRef();
    Foam::Field<TVector> & Uf_ = Uf.primitiveFieldRef();
    Foam::Field<TScalar> & Pf_ = Pf.primitiveFieldRef();
    Foam::Field<TScalar> & Ef_ = Ef.primitiveFieldRef();
    Foam::Field<TScalar> & Hthermodynamicalf_ = Hthermodynamicalf.primitiveFieldRef();
    Foam::Field<TScalar> & MolWeightf_ = MolWeightf.primitiveFieldRef();
    Foam::Field<TScalar> & Inductionf_ = Inductionf.primitiveFieldRef();
    const Foam::GeometricField<TVector, fvsPatchField, surfaceMesh> & _Sf = mesh.Sf();
    const Foam::GeometricField<TScalar, fvsPatchField, surfaceMesh> & _magSf = mesh.magSf();
    const Foam::GeometricField<TVector, fvPatchField, volMesh> & CellCentres = mesh.C();
    const Foam::GeometricField<TVector, fvsPatchField, surfaceMesh> & FaceCentres = mesh.Cf();

    for (Foam::label fi = 0; fi < P.size(); fi++)
    {
        TVector Normal = _Sf[fi] / _magSf[fi];
		// Ячейка №1 - P
        TScalar dP = mag(FaceCentres[fi] - CellCentres[P[fi]]);
        TScalar rho1 = Rho_[P[fi]];
		TScalar p1 = p_[P[fi]];
        TVector u1 = U_[P[fi]];
        TScalar MolWeight1 = MolWeight_[P[fi]];
        TScalar Induction1 = Induction_[P[fi]];
        TScalar c1 = c_[P[fi]];
        TScalar gamma1 = gamma_[P[fi]];
        TScalar T1 = T_[P[fi]];
        //----------------------------------------------------------------
        // Разложение скорости
        TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		TVector UT1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
        
		// Ячейка №2 - N
        TScalar dN =  mag(CellCentres[N[fi]] - FaceCentres[fi]);
        TScalar rho2 = Rho_[N[fi]];
		TScalar p2 = p_[N[fi]];
        TVector u2 = U_[N[fi]];
        TScalar MolWeight2 = MolWeight_[N[fi]];
        TScalar Induction2 = Induction_[N[fi]];
        TScalar c2 = c_[N[fi]];
        TScalar gamma2 = gamma_[N[fi]];
        TScalar T2 = T_[N[fi]];
        //----------------------------------------------------------------
		// Разложение скорости
		TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
		TVector UT2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
        
        //----------------------------------------------------------------
		// out
		TScalar PP;
		TScalar RR;
		TScalar UU;
		TScalar Ustar;                      // скорость контактного разрыва
		TScalar W = 0.0;
		int sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
        if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
		Rhof_[fi] = RR;
        Pf_[fi] = PP;
		if(sln < 0)
		{
            Uf_[fi] = UT1 + (Normal * UU);
        	MolWeightf_[fi] = MolWeight1;
        	Inductionf_[fi] = Induction1;
            ChemistryElement Element1(Induction1, MolWeight1);
            Hthermodynamicalf_[fi] = Thermo.Hthermodynamical(PP, RR, Element1);
            Ef_[fi] = Thermo.Hpotential(PP, RR, Element1);
            continue;
		}
        if(sln > 0)
		{
            Uf_[fi] = UT2 + (Normal * UU);
        	MolWeightf_[fi] = MolWeight2;
        	Inductionf_[fi] = Induction2;
            ChemistryElement Element2(Induction2, MolWeight2);
            Hthermodynamicalf_[fi] = Thermo.Hthermodynamical(PP, RR, Element2);
            Ef_[fi] = Thermo.Hpotential(PP, RR, Element2);
            continue;
		}
        // Ошибка: sln=0
        Info << "ERROR!" << endl;
        FatalError.abort();
        //Error
    }

    // Interpolate across coupled patches using given lamb  das
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Rhobf = Rhof.boundaryFieldRef();
    typename Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Ubf = Uf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Pbf = Pf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Ebf = Ef.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Hthermodynamicalbf = Hthermodynamicalf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & MolWeightbf = MolWeightf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Inductionbf = Inductionf.boundaryFieldRef();

    forAll(mesh.boundaryMesh(), ipatch)
    {
        Foam::word BCtype = mesh.boundaryMesh().types()[ipatch];
        const Foam::fvsPatchScalarField & pLambda = lambdas.boundaryField()[ipatch];
        const typename SFType::Patch & pSf = Sf.boundaryField()[ipatch];

        Foam::fvsPatchField<TScalar> & Rhosf = Rhobf[ipatch];
        Foam::fvsPatchField<TVector> & Usf = Ubf[ipatch];
        Foam::fvsPatchField<TScalar> & Psf = Pbf[ipatch];
        Foam::fvsPatchField<TScalar> & Esf = Ebf[ipatch];
        Foam::fvsPatchField<TScalar> & Hthermodynamicalsf = Hthermodynamicalbf[ipatch];
        Foam::fvsPatchField<TScalar> & MolWeightsf = MolWeightbf[ipatch];
        Foam::fvsPatchField<TScalar> & Inductionsf = Inductionbf[ipatch];

        // 2 в 1:
        // 1 - межпроцессорное взаимодействие processor
        // 2 - периодическое граничное условие processorPatch
        bool UseProcessorBounary = ((TimeWall && BCtype == "processor") || (!TimeWall && isA<processorFvPatch>(mesh.boundary()[ipatch])));
        // граничное условия для двух областей, но обе области нахидятся под одним процессором
        bool UsePatchBounary = (BCtype != "processor") && (BCtype != "cyclic") && (BCtype != "processorCyclic");
        // Если межпроцессорное взаимодействие
        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && UseProcessorBounary)
        {
            const Foam::fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const Foam::fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const Foam::fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const Foam::fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча P - Internal - Полоса ячеек этого процессора
            Foam::tmp<Foam::Field<TScalar>> tmpRhoInternal(Rho.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TVector>> tmpUInternal(U.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpPInternal(p.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpMolWeightInternal(MolWeight.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpInductionInternal(Induction.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpcInternal(c.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpgammaInternal(gamma.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpEInternal(E.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpTInternal(T.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpdeltaInternal = mag(FaceCentres_ - CellCentres_);
            // Внутренние ссылки
            Foam::Field<TScalar> & RhoInternal = tmpRhoInternal.ref();
            Foam::Field<TVector> & UInternal = tmpUInternal.ref();
            Foam::Field<TScalar> & PInternal = tmpPInternal.ref();
            Foam::Field<TScalar> & EInternal = tmpEInternal.ref();
            Foam::Field<TScalar> & TInternal = tmpTInternal.ref();
            Foam::Field<TScalar> & MolWeightInternal = tmpMolWeightInternal.ref();
            Foam::Field<TScalar> & InductionInternal = tmpInductionInternal.ref();
            Foam::Field<TScalar> & gammaInternal = tmpgammaInternal.ref();
            Foam::Field<TScalar> & cInternal = tmpcInternal.ref();
            Foam::Field<TScalar> & deltaInternal = tmpdeltaInternal.ref();

            // Для конкретного патча N - Neighbour - Полоса ячеек соседнего процессора
            tmp<Field<TScalar>> tmpRhoNeighbour(Rho.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TVector>> tmpUNeighbour(U.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpPNeighbour(p.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpENeighbour(E.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpTNeighbour(T.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpMolWeightNeighbour(MolWeight.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpInductionNeighbour(Induction.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpgammaNeighbour(gamma.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpcNeighbour(c.boundaryField()[ipatch].patchInternalField());
            tmp<Field<TScalar>> tmpdeltaNeighbour = mag(FaceCentres_ - CellCentres_);
            // Внутренние ссылки
            Field<TScalar> & RhoNeighbour = tmpRhoNeighbour.ref();
            Field<TVector> & UNeighbour = tmpUNeighbour.ref();
            Field<TScalar> & PNeighbour = tmpPNeighbour.ref();
            Field<TScalar> & ENeighbour = tmpENeighbour.ref();
            Field<TScalar> & TNeighbour = tmpTNeighbour.ref();
            Field<TScalar> & MolWeightNeighbour = tmpMolWeightNeighbour.ref();
            Field<TScalar> & InductionNeighbour = tmpInductionNeighbour.ref();
            Field<TScalar> & gammaNeighbour = tmpgammaNeighbour.ref();
            Field<TScalar> & cNeighbour = tmpcNeighbour.ref();
            Field<TScalar> & deltaNeighbour = tmpdeltaNeighbour.ref();
            
            const Foam::processorPolyPatch & ProcessorPolyPatch = refCast<const Foam::processorPolyPatch>(mesh.boundaryMesh()[ipatch]);
            MPIexchange<TScalar, TVector>
            (
                RhoInternal, UInternal, PInternal, EInternal, TInternal, MolWeightInternal, InductionInternal, gammaInternal, cInternal, deltaInternal,
                RhoNeighbour, UNeighbour, PNeighbour, ENeighbour, TNeighbour,  MolWeightNeighbour, InductionNeighbour, gammaNeighbour, cNeighbour, deltaNeighbour,
                ProcessorPolyPatch
            );
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar dP = deltaInternal[i];
                TScalar rho1 = RhoInternal[i];
                TScalar p1 = PInternal[i];
                TVector u1 = UInternal[i];
                TScalar MolWeight1 = MolWeightInternal[i];
                TScalar Induction1 = InductionInternal[i];
                TScalar c1 = cInternal[i];
                TScalar gamma1 = gammaInternal[i];
                TScalar T1 = TInternal[i];
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector UT1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                // Ячейка №2 - N
                TScalar dN = deltaNeighbour[i];
                TScalar rho2 = RhoNeighbour[i];
                TScalar p2 = PNeighbour[i];
                TVector u2 = UNeighbour[i];
                TScalar MolWeight2 = MolWeightNeighbour[i];
                TScalar Induction2 = InductionNeighbour[i];
                TScalar c2 = cNeighbour[i];
                TScalar gamma2 = gammaNeighbour[i];
                TScalar T2 = TNeighbour[i];
                // Поправка Колгана для ячейки №2 - N
                // Разложение скорости
                TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
                TVector UT2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
                //----------------------------------------------------
                // out
                TScalar PP;
                TScalar RR;
                TScalar UU;
                TScalar Ustar;                      // скорость контактного разрыва
                TScalar W = 0.0;
                int sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
                if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
                Rhosf[i] = RR;
                Psf[i] = PP;
                if(sln < 0)
                {
                    Usf[i] = UT1 + (Normal * UU);
                    MolWeightsf[i] = MolWeight1;
                    Inductionsf[i] = Induction1;
                    ChemistryElement Element1(Induction1, MolWeight1);
                    Hthermodynamicalsf[i] = Thermo.Hthermodynamical(PP, RR, Element1);
                    Esf[i] = Thermo.Hpotential(PP, RR, Element1);
                    continue;
                }
                if(sln > 0)
                {
                    Usf[i] = UT2 + (Normal * UU);
                    MolWeightsf[i] = MolWeight2;
                    Inductionsf[i] = Induction2;
                    ChemistryElement Element2(Induction2, MolWeight2);
                    Hthermodynamicalsf[i] = Thermo.Hthermodynamical(PP, RR, Element2);
                    Esf[i] = Thermo.Hpotential(PP, RR, Element2);
                    continue;
                }
                // Ошибка: sln=0
                Info << "ERROR!" << endl;
                FatalError.abort();
                //Error
            }
        }else if (mesh.boundaryMesh()[ipatch].coupled() && UsePatchBounary) // PatchIndex
        {
            const Foam::fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const Foam::fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const Foam::fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const Foam::fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча поля для P
            Foam::tmp<Foam::Field<TScalar>> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TVector>> UInternal = U.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> PInternal = p.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> cInternal = c.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> TInternal = T.boundaryField()[ipatch].patchInternalField();

            // Для конкретного патча параметры N
            Foam::tmp<Foam::Field<TScalar>> RhoNeighbour = Rho.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TVector>> UNeighbour = U.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> PNeighbour = p.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> MolWeightNeighbour = MolWeight.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> InductionNeighbour = Induction.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> cNeighbour = c.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> gammaNeighbour = gamma.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> TNeighbour = T.boundaryField()[ipatch].patchNeighbourField();
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar dP = mag(FaceCentres_[i] - CellCentres_[i]);
                TScalar rho1 = RhoInternal.ref()[i];
                TScalar p1 = PInternal.ref()[i];
                TVector u1 = UInternal.ref()[i];
                TScalar MolWeight1 = MolWeightInternal.ref()[i];
                TScalar Induction1 = InductionInternal.ref()[i];
                TScalar c1 = cInternal.ref()[i];
                TScalar gamma1 = gammaInternal.ref()[i];
                TScalar T1 = TInternal.ref()[i];
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector UT1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                // Ячейка №2 - N
                TScalar dN = mag(CellCentres_[i] - FaceCentres_[i]);
                TScalar rho2 = RhoNeighbour.ref()[i];
                TScalar p2 = PNeighbour.ref()[i];
                TVector u2 = UNeighbour.ref()[i];
                TScalar MolWeight2 = MolWeightNeighbour.ref()[i];
                TScalar Induction2 = InductionNeighbour.ref()[i];
                TScalar c2 = cNeighbour.ref()[i];
                TScalar gamma2 = gammaNeighbour.ref()[i];
                TScalar T2 = TNeighbour.ref()[i];
                // Поправка Колгана для ячейки №2 - N
                // Разложение скорости
                TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
                TVector UT2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
                //----------------------------------------------------
                // out
                TScalar PP;
                TScalar RR;
                TScalar UU;
                TScalar Ustar;                      // скорость контактного разрыва
                TScalar W = 0.0;
                int sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
                if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
                Rhosf[i] = RR;
                Psf[i] = PP;
                if(sln < 0)
                {
                    Usf[i] = UT1 + (Normal * UU);
                    MolWeightsf[i] = MolWeight1;
                    Inductionsf[i] = Induction1;
                    ChemistryElement Element1(Induction1, MolWeight1);
                    Hthermodynamicalsf[i] = Thermo.Hthermodynamical(PP, RR, Element1);
                    Esf[i] = Thermo.Hpotential(PP, RR, Element1);
                    continue;
                }
                if(sln > 0)
                {
                    Usf[i] = UT2 + (Normal * UU);
                    MolWeightsf[i] = MolWeight2;
                    Inductionsf[i] = Induction2;
                    ChemistryElement Element2(Induction2, MolWeight2);
                    Hthermodynamicalsf[i] = Thermo.Hthermodynamical(PP, RR, Element2);
                    Esf[i] = Thermo.Hpotential(PP, RR, Element2);
                    continue;
                }
                // Ошибка: sln=0
                Info << "ERROR!" << endl;
                FatalError.abort();
                //Error
            }
        }
        else
        {
            forAll(Rho.boundaryField()[ipatch], i)
            {
                Rhosf[i] = Rho.boundaryField()[ipatch][i];
                Usf[i] = U.boundaryField()[ipatch][i];
                Psf[i] = p.boundaryField()[ipatch][i];
                Esf[i] = E.boundaryField()[ipatch][i];
                MolWeightsf[i] = MolWeight.boundaryField()[ipatch][i];
                Inductionsf[i] = Induction.boundaryField()[ipatch][i];
                ChemistryElement Element(Inductionsf[i], MolWeightsf[i]);
                Hthermodynamicalsf[i] = Thermo.Hthermodynamical(Psf[i], Rhosf[i], Element);
            }
        }
    }
    tlambdas.clear();
}
