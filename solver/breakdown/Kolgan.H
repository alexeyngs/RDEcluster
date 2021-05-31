template<typename TScalar, typename TYPE>
inline TYPE gradient(const TScalar dH, const TYPE ff1, const TYPE ff2, const TYPE ff3)
// dH - расстояние от центра до грани - половина размера ячейки
{
    TYPE outresult = (ff3 - ff1)/(dH*2.0);
    if(Foam::mag(outresult) > Foam::mag(ff3 - ff2)/dH) outresult = (ff3 - ff2)/dH;
    if(Foam::mag(outresult) > Foam::mag(ff2 - ff1)/dH) outresult = (ff2 - ff1)/dH;
    return outresult;
};

template<typename TScalar, typename TYPE>
inline TYPE gradient(const TScalar L1, const TScalar L2, const TYPE ff1, const TYPE ff2, const TYPE ff3)
// L1 - размер левой ячейки
// L2 - размер правой ячейки
{
    TYPE outresult = (ff3 - ff1)/((L1 + L2) / 2.0);
    if(Foam::mag(outresult) > Foam::mag(ff3 - ff2) / (L2 / 2.0)) outresult = (ff3 - ff2) / (L2 / 2.0);
    if(Foam::mag(outresult) > Foam::mag(ff2 - ff1) / (L1 / 2.0)) outresult = (ff2 - ff1) / (L1 / 2.0);
    return outresult;
};

inline Foam::label OpposideFaceLabel(const fvMesh & mesh, Foam::label celli, Foam::label face)
{
    // Координата текущей грани
    const point & facevector = mesh.faceCentres()[celli];
    // Список всех граней
    const Foam::labelList & cFaces = mesh.cells()[celli];
    Foam::label result = 0;
    Foam::scalar maxdelta = 0.0;
    forAll(cFaces, cFacei)
    {
        Foam::label facei = cFaces[cFacei]; // получаем индекс сетки
        const Foam::point & p = mesh.faceCentres()[facei];  // центр ячейки с определеннаи индексом
        Foam::scalar delta = Foam::mag(p - facevector);
        if(delta > maxdelta)
        {
            delta = maxdelta;
            result = facei;
        }
    }
    return result;
};

template<typename TYPE>
inline TYPE Face(const Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh> & Field, Foam::label index)
{
    if(index < Field.size()) return Field[index];
    else
    {
        index -= Field.size();
        forAll(Field.boundaryField(), patchi)
        {
            if(index >= Field.boundaryField()[patchi].size())
            {
                index -= Field.boundaryField()[patchi].size();
            }else
            {
                return Field.boundaryField()[patchi][index];
            }
        }
    }
    return Field[0];
};
//=================================================================================================================
template<typename TScalar, typename TVector, typename SFType>
inline void Kolgan
(
    rhoRDEThermo & Thermo,
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
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rholf,
    const Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Ulf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Plf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Elf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightlf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionlf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rhof,
    Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Uf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Pf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ef,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Hthermodynamicalf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionf
)
{
    const fvMesh & mesh = Rho.mesh();
    const labelUList & P = mesh.owner();
    const labelUList & N = mesh.neighbour();
    //-------------------------------------------------------------------------------
    const tmp<surfaceScalarField> & tlambdas = Rho.mesh().surfaceInterpolation::weights();
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
        const TVector Normal = _Sf[fi] / _magSf[fi];
		// Ячейка №1 - P
        TScalar dP = mag(FaceCentres[fi] - CellCentres[P[fi]]);
        TScalar rho1 = Rho_[P[fi]];
		TScalar p1 = p_[P[fi]];
        TVector u1 = U_[P[fi]];
        TScalar E1 = E_[P[fi]];
        TScalar MolWeight1 = MolWeight_[P[fi]];
        TScalar Induction1 = Induction_[P[fi]];
        TScalar c1 = c_[P[fi]];
        TScalar gamma1 = gamma_[P[fi]];
        TScalar T1 = T_[P[fi]];
        // Поправка Колгана для ячейки №1 - P
        Foam::label OpposideP = OpposideFaceLabel(mesh, P[fi], fi); // Индекс противоположной грани
        rho1 += gradient(dP, Face(Rholf, OpposideP), rho1, Rholf[fi])*dP;
        p1 += gradient(dP, Face(Plf, OpposideP), p1, Plf[fi])*dP;
        u1 += gradient(dP, Face(Ulf, OpposideP), u1, Ulf[fi])*dP;
        E1 += gradient(dP, Face(Elf, OpposideP), E1, Elf[fi])*dP;
        MolWeight1 += gradient(dP, Face(MolWeightlf, OpposideP), MolWeight1, MolWeightlf[fi])*dP;
        Induction1 += gradient(dP, Face(Inductionlf, OpposideP), Induction1, Inductionlf[fi])*dP;
        ChemistryElement Element1(Induction1, MolWeight1);
        // Разложение скорости
        TScalar U1 = u1 & Normal;           // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		TVector UT1 = u1 - (Normal * U1);    // Тангенсальная компонента скорости в ячейке 1

        //----------------------------------------------------------------
		// Ячейка №2 - N
        TScalar dN =  mag(CellCentres[N[fi]] - FaceCentres[fi]);
        TScalar rho2 = Rho_[N[fi]];
		TScalar p2 = p_[N[fi]];
        TVector u2 = U_[N[fi]];
        TScalar E2 = E_[N[fi]];
        TScalar MolWeight2 = MolWeight_[N[fi]];
        TScalar Induction2 = Induction_[N[fi]];
        TScalar c2 = c_[N[fi]];
        TScalar gamma2 = gamma_[N[fi]];
        TScalar T2 = T_[N[fi]];
        // Поправка Колгана для ячейки №2 - N
        Foam::label OpposideN = OpposideFaceLabel(mesh, N[fi], fi); // Индекс противоположной грани
        rho2 -= gradient(dN, Rholf[fi], rho2, Face(Rholf, OpposideN))*dN;
        p2 -= gradient(dN, Plf[fi], p2, Face(Plf, OpposideN))*dN;
        u2 -= gradient(dN, Ulf[fi], u2, Face(Ulf, OpposideN))*dN;
        E2 -= gradient(dN, Elf[fi], E2, Face(Elf, OpposideN))*dN;
        MolWeight2 -= gradient(dN, MolWeightlf[fi], MolWeight2, Face(MolWeightlf, OpposideN))*dN;
        Induction2 -= gradient(dN, Inductionlf[fi], Induction2, Face(Inductionlf, OpposideN))*dN;
        ChemistryElement Element2(Induction2, MolWeight2);
		// Разложение скорости
		TScalar U2 = u2 & Normal;           // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
		TVector UT2 = u2 - (Normal * U2);    // Тангенсальная компонента скорости в ячейке 2
        
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
        Info << "ERROR KOLGAN! TRY GODUNOV" << endl;
        rho1 = Rho_[P[fi]];
        p1 = p_[P[fi]];
        u1 = U_[P[fi]];
        E1 = E_[P[fi]];
        MolWeight1 = MolWeight_[P[fi]];
        Induction1 = Induction_[P[fi]];
        U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
        UT1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
        c1 = c_[P[fi]];
        gamma1 = gamma_[P[fi]];
        T1 = T_[P[fi]];
        //---------------------------------
        rho2 = Rho_[N[fi]];
        p2 = p_[N[fi]];
        u2 = U_[N[fi]];
        E2 = E_[N[fi]];
        MolWeight2 = MolWeight_[N[fi]];
        Induction2 = Induction_[N[fi]];
        U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
        UT2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 1
        c2 = c_[N[fi]];
        gamma2 = gamma_[N[fi]];
        T2 = T_[N[fi]];
        //----------------------------------
        sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
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
    // Interpolate across coupled patches using given lambdas
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Rhobf = Rhof.boundaryFieldRef();
    typename GeometricField<TVector, fvsPatchField, surfaceMesh>::Boundary & Ubf = Uf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Pbf = Pf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Ebf = Ef.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Hthermodynamicalbf = Hthermodynamicalf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & MolWeightbf = MolWeightf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Inductionbf = Inductionf.boundaryFieldRef();

    forAll(mesh.boundaryMesh(), ipatch)
    {
        word BCtype = mesh.boundaryMesh().types()[ipatch];
        const fvsPatchScalarField & pLambda = lambdas.boundaryField()[ipatch];
        const typename SFType::Patch & pSf = Sf.boundaryField()[ipatch];

        fvsPatchField<TScalar> & Rhosf = Rhobf[ipatch];
        fvsPatchField<TVector> & Usf = Ubf[ipatch];
        fvsPatchField<TScalar> & Psf = Pbf[ipatch];
        fvsPatchField<TScalar> & Esf = Ebf[ipatch];
        fvsPatchField<TScalar> & Hthermodynamicalsf = Hthermodynamicalbf[ipatch];
        fvsPatchField<TScalar> & MolWeightsf = MolWeightbf[ipatch];
        fvsPatchField<TScalar> & Inductionsf = Inductionbf[ipatch];

        // 2 в 1:
        // 1 - межпроцессорное взаимодействие processor
        // 2 - периодическое граничное условие processorPatch
        bool UseProcessorBounary = ((TimeWall && BCtype == "processor") || (!TimeWall && isA<processorFvPatch>(mesh.boundary()[ipatch])));
        // граничное условия для двух областей, но обе области нахидятся под одним процессором
        bool UsePatchBounary = (BCtype != "processor") && (BCtype != "cyclic") && (BCtype != "processorCyclic");
        // Если межпроцессорное взаимодействие
        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && UseProcessorBounary)
        {
            const fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];
            // Список индексов граней на границе
            const Foam::labelList & boundaryFaces = mesh.boundaryMesh()[ipatch].meshEdges();
            // delta "свой"
            Foam::tmp<Foam::Field<TScalar>> tmpdeltaInternal = mag(FaceCentres_ - CellCentres_);
            Foam::Field<TScalar> & deltaInternal = tmpdeltaInternal.ref();
            // delta "соседнего процессора"
            Foam::tmp<Foam::Field<TScalar>> tmpdeltaNeighbour = mag(FaceCentres_ - CellCentres_);
            Foam::Field<TScalar> & deltaNeighbour = tmpdeltaNeighbour.ref();
            // Список индексов ячеек на границе(пограничные ячейки)
            const labelList & InternalCells = mesh.C().boundaryField()[ipatch].patch().patch().faceCells();
            // Получение идентификатора процессора(себя и соседнего)
            const Foam::processorPolyPatch & ProcessorPolyPatch = refCast<const Foam::processorPolyPatch>(mesh.boundaryMesh()[ipatch]);
            //--------------------------------------------------------------------------------------------------------
            // Для конкретного патча P - Internal - Полоса ячеек этого процессора
            Foam::tmp<Foam::Field<TScalar>> tmpRhoInternal(Rho.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TVector>> tmpUInternal(U.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpPInternal(p.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpEInternal(E.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpTInternal(T.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpMolWeightInternal(MolWeight.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpInductionInternal(Induction.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpgammaInternal(gamma.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpcInternal(c.boundaryField()[ipatch].patchInternalField());
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
            //--------------------------------------------------------------------------------------------------------
            // Для конкретного патча N - экспортируем полоски ячеек в соседний процессор
            Foam::tmp<Foam::Field<TScalar>> tmpRhoExport(Rho.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TVector>> tmpUExport(U.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpPExport(p.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpEExport(E.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpTExport(T.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpMolWeightExport(MolWeight.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpInductionExport(Induction.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpgammaExport(gamma.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpcExport(c.boundaryField()[ipatch].patchInternalField());
            // Внутренние ссылки
            Foam::Field<TScalar> & RhoExport = tmpRhoExport.ref();
            Foam::Field<TVector> & UExport = tmpUExport.ref();
            Foam::Field<TScalar> & PExport = tmpPExport.ref();
            Foam::Field<TScalar> & EExport = tmpEExport.ref();
            Foam::Field<TScalar> & TExport = tmpTExport.ref();
            Foam::Field<TScalar> & MolWeightExport = tmpMolWeightExport.ref();
            Foam::Field<TScalar> & InductionExport = tmpInductionExport.ref();
            Foam::Field<TScalar> & gammaExport = tmpgammaExport.ref();
            Foam::Field<TScalar> & cExport = tmpcExport.ref();
            //--------------------------------------------------------------------------------------------------------
            // Для конкретного патча N - ипортируем полоски ячеек из соседнего процессора     
            Foam::tmp<Foam::Field<TScalar>> tmpRhoImport(Rho.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TVector>> tmpUImport(U.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpPImport(p.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpEImport(E.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpTImport(T.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpMolWeightImport(MolWeight.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpInductionImport(Induction.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpgammaImport(gamma.boundaryField()[ipatch].patchInternalField());
            Foam::tmp<Foam::Field<TScalar>> tmpcImport(c.boundaryField()[ipatch].patchInternalField());
            // Внутренние ссылки
            Foam::Field<TScalar> & RhoImport = tmpRhoImport.ref();
            Foam::Field<TVector> & UImport = tmpUImport.ref();
            Foam::Field<TScalar> & PImport = tmpPImport.ref();
            Foam::Field<TScalar> & EImport = tmpEImport.ref();
            Foam::Field<TScalar> & TImport = tmpTImport.ref();
            Foam::Field<TScalar> & MolWeightImport = tmpMolWeightImport.ref();
            Foam::Field<TScalar> & InductionImport = tmpInductionImport.ref();
            Foam::Field<TScalar> & gammaImport = tmpgammaImport.ref();
            Foam::Field<TScalar> & cImport = tmpcImport.ref();
            //--------------------------------------------------------------------------------------------------------
            // Полоски(патчи) грани от линейной интерполяции(lf) - slf
            const fvsPatchField<TScalar> & Rhoslf = Rholf.boundaryField()[ipatch];
            const fvsPatchField<TVector> & Uslf = Ulf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Pslf = Plf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Eslf = Elf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & MolWeightslf = MolWeightlf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Inductionslf = Inductionlf.boundaryField()[ipatch];
            //--------------------------------------------------------------------------------------------------------
            // Заполнение поправок
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TScalar delta = deltaInternal[i];
                // Индекс пограничной грани
                Foam::label Edge = boundaryFaces[i];
                // Индекс противоположной грани(от пограничной)
                Foam::label Opposide = OpposideFaceLabel(mesh, InternalCells[i], Edge);
                // Собственник - свои ячейки
                RhoInternal[i] += gradient(delta, Face(Rholf, Opposide), RhoInternal[i], Rhoslf[i])*delta;
                UInternal[i] += gradient(delta, Face(Ulf, Opposide), UInternal[i], Uslf[i])*delta;
                PInternal[i] += gradient(delta, Face(Plf, Opposide), PInternal[i], Pslf[i])*delta;
                EInternal[i] += gradient(delta, Face(Elf, Opposide), EInternal[i], Eslf[i])*delta;
                MolWeightInternal[i] += gradient(delta, Face(MolWeightlf, Opposide), MolWeightInternal[i], MolWeightslf[i])*delta;
                InductionInternal[i] += gradient(delta, Face(Inductionlf, Opposide), InductionInternal[i], Inductionslf[i])*delta;
                // Соседние ячейка на экспорт
                RhoExport[i] -= gradient(delta, Rhoslf[i], RhoInternal[i], Face(Rholf, Opposide))*delta;
                UExport[i] -= gradient(delta, Uslf[i], UInternal[i], Face(Ulf, Opposide))*delta;
                PExport[i] -= gradient(delta, Pslf[i], PInternal[i], Face(Plf, Opposide))*delta;
                EExport[i] -= gradient(delta, Eslf[i], EInternal[i], Face(Elf, Opposide))*delta;
                MolWeightExport[i] -= gradient(delta, MolWeightslf[i], MolWeightInternal[i], Face(MolWeightlf, Opposide))*delta;
                InductionExport[i] -= gradient(delta, Inductionslf[i], InductionInternal[i], Face(Inductionlf, Opposide))*delta;
                // Проверка на отрицательность строго-положительных параметров
                if(RhoInternal[i] < 0.1) RhoInternal[i] = Rho.boundaryField()[ipatch].patchInternalField().ref()[i];
                if(RhoExport[i] < 0.1) RhoExport[i] = Rho.boundaryField()[ipatch].patchInternalField().ref()[i];
                if(PInternal[i] < 1.0) PInternal[i] = p.boundaryField()[ipatch].patchInternalField().ref()[i];
                if(PExport[i] < 1.0) PExport[i] = p.boundaryField()[ipatch].patchInternalField().ref()[i];
            }
            // Обмен данных с соседним процессором(отдаем свои - Export и получаем от соседа - Import)
            MPIexchange<TScalar, TVector>
            (
                RhoExport, UExport, PExport, EExport, TExport, MolWeightExport, InductionExport, gammaExport, cExport, deltaInternal,
                RhoImport, UImport, PImport, EImport, TImport, MolWeightImport, InductionImport, gammaImport, cImport, deltaNeighbour,
                ProcessorPolyPatch
            );
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar rho1 = RhoInternal[i];
                TScalar p1 = PInternal[i];
                TVector u1 = UInternal[i];
                TScalar gamma1 = gammaInternal[i];
                TScalar E1 = EInternal[i];
                TScalar T1 = TInternal[i];
                TScalar MolWeight1 = MolWeightInternal[i];
                TScalar Induction1 = InductionInternal[i];
                ChemistryElement Element1(Induction1, MolWeight1);
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector UT1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                TScalar c1 = cInternal[i];
                // Ячейка №2 - N
                TScalar rho2 = RhoImport[i];
                TScalar p2 = PImport[i];
                TVector u2 = UImport[i];
                TScalar gamma2 = gammaImport[i];
                TScalar E2 = EImport[i];
                TScalar T2 = TImport[i];
                TScalar MolWeight2 = MolWeightImport[i];
                TScalar Induction2 = InductionImport[i];
                ChemistryElement Element2(Induction2, MolWeight2);
                // Поправка Колгана для ячейки №2 - N
                // Разложение скорости
                TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
                TVector UT2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
                TScalar c2 = cImport[i];
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
                Foam::Info << "ERROR in Kolgan MPI exchange!" << Foam::endl;
                Foam::Info << "sln = 0: p1 = " << p1 << ", pho1 = " << rho1 << ", U1 = " << U1 << ", c1 = " << c1 << Foam::endl;
                Foam::Info << "sln = 0: p2 = " << p2 << ", pho2 = " << rho2 << ", U2 = " << U2 << ", c2 = " << c2 << Foam::endl;
                Foam::Info << "sln = 0: PP = " << PP << ", RR = " << RR << ", UU = " << UU << Foam::endl;
                FatalError.abort();
                //Error
            }
        }else if (mesh.boundaryMesh()[ipatch].coupled() && UsePatchBounary) // PatchIndex
        {
            goto patch;
        }
        else
        {
            patch:
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
