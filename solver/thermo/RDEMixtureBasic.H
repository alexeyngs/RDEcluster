#ifndef RDEMixtureBasic_H
#define RDEMixtureBasic_H

//----------------------------------------------------------------
const Foam::scalar Calorie = 4.1868;
const Foam::scalar KiloCalorie = 4186.8;
const Foam::scalar Kplus = 6.0e8;                       // K+
const Foam::scalar T0 = 300.0;                          // T0
const Foam::scalar Tinert = 300.0;                      // Температура, ниже которой хим. кинетика запрещена

template<class ThermoType>
inline Foam::scalar GetMoles(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::word & name)
{
    if(!species.found(name)) return 0.0;
    if(!speciesData(species[name])) return 0.0;
    return speciesData[species[name]].moles();
};

template<class ThermoType>
inline Foam::scalar GetEdissociation(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::word & name)
{
    if(!species.found(name)) return 0.0;
    if(!speciesData(species[name])) return 0.0;
    return speciesData[species[name]].Edissociation();
};

template<class ThermoType>
inline Foam::scalar GetTheta(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::word & name)
{
    if(!species.found(name)) return 0.0;
    if(!speciesData(species[name])) return 0.0;
    return speciesData[species[name]].theta();
};

//----------------------------------------------------------------
namespace Foam{
template<class ThermoType>
class RDEMixtureBasic : public IMixture, public speciesTable, public multiComponentMixture<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureBasic(const RDEMixtureBasic &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureBasic &) = delete;
public:
    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;
    //- Table of species composition
    speciesCompositionTable speciesComposition_;
    volScalarField MolWeight_;
    volScalarField FMolWeight_;
    volScalarField Induction_;
    volScalarField FInduction_;
    scalar powerN;
    scalar theta;
    scalar beta;
    scalar Einduction;          // Энергия активации [Дж/моль]
    scalar Kinduction;          // Константа активации в моль*с/м3
    scalar StartInduction;
    scalar CompareWeight;
    // CJ
    scalar CJ_D;
    scalar CJ_P;
    scalar CJ_T;
    scalar CJ_Rho;
    scalar CJ_MolWeight;
    scalar CJ_Kminus;
    //----------------------------------------------------------------
    scalar MolWeight0;
    scalar MolWeightA;          // Атомарная масса
    scalar MolWeightMin;        // Масса предельной диссоциации
    scalar MolWeightMax;        // Масса предельной рекомбинации
    scalar Edissociation;       // Энергия диссоциации - коэфициент
    scalar Qmax;                // Энергия
    scalar fuel;                // Количество атомов кислорода для полного сгорания
    scalar Sigma;               // Суммарная масса всей смеси
    scalar atoms;               // Общее количество атомов с учетом молей
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureBasic(const dictionary &, const fvMesh &, const word &);
    //- Destructor
    virtual ~RDEMixtureBasic();
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureBasic<" + ThermoType::typeName() + '>';
    }
    //- Read dictionary
    void read(const dictionary&);
    //- Table of species composition
    const speciesCompositionTable & specieComposition() const
    {
        return speciesComposition_;
    }
    Foam::scalar W(const label speciei) const
    {
        const ThermoType & type = multiComponentMixture<ThermoType>::getLocalThermo(speciei);
        return type.W()/type.Y();
    }
    void Inert(const label celli);
    void CalculateCJ();
    scalar GetMoles(const word & name)const;
    const word InfoElement(const label celli);
    const word InfoElement(const label patchi, const label facei);
    //--------------------------------------------------------------------------------------------------
    virtual scalar GetA(const scalar MolWeight, const scalar T) const override;
    virtual scalar Cv(const scalar P, const scalar T, const scalar MolWeight)const override;
    virtual scalar GetTmax(const Foam::scalar MolWeight, const scalar Uthermodinamical) const override;
};
} // End namespace Foam
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------

template<class ThermoType>
Foam::RDEMixtureBasic<ThermoType>::RDEMixtureBasic(const Foam::dictionary & thermoDict, const Foam::fvMesh & mesh, const Foam::word & phaseName):
    speciesTable(),
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    speciesComposition_(),
    MolWeight_
    (
        IOobject
        (
            "MolWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    FMolWeight_
    (
        IOobject
        (
            "FMolWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("dimVelocity*dimMoles", dimensionSet(1, -3, -1, 0, 1, 0, 0), 0.0)
    ),
    Induction_
    (
        IOobject
        (
            "Induction",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    FInduction_
    (
        IOobject
        (
            "FInduction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("dimVelocity", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    )
{
    const Foam::Time & runTime = mesh.time();
    const ThermoType & mixture = multiComponentMixture<ThermoType>::cellMixture(0);
    //-----------------------------------------------------------------
    MolWeight0 = mixture.W();
    MolWeightA = mixture.W()/mixture.atoms();           // Общее количество атомов с учетом молей
    MolWeightMin = mixture.W()/mixture.decomposition(); // Масса предельной диссоциации
    fuel = mixture.needO()*mixture.moles();             // Количество атомов кислорода для полного сгорания
    Sigma = mixture.W()*mixture.moles();                // Суммарная масса всей смеси
    atoms = mixture.atoms()*mixture.moles();            // Общее количество атомов с учетом молей
    // Получение параметров, которые зависят от состава смеси
    powerN = readScalar(thermoDict.subDict("thermo").lookup("powerN"));
    beta = readScalar(thermoDict.subDict("thermo").lookup("beta"));
    StartInduction = readScalar(thermoDict.subDict("thermo").lookup("StartInduction"));
    theta = thermoDict.subDict("thermo").lookupOrDefault("theta", 0.0);
    Qmax = thermoDict.subDict("thermo").lookupOrDefault("Qmax", 0.0);
    CompareWeight = thermoDict.subDict("thermo").lookupOrDefault("CompareWeight", -1.0);
    // Параметры Чепмена-Жуге
    CJ_D = readScalar(thermoDict.subDict("CJ").lookup("D"));
    CJ_P = readScalar(thermoDict.subDict("CJ").lookup("P"));
    CJ_T = readScalar(thermoDict.subDict("CJ").lookup("T"));
    //-----------------
    CJ_Rho = thermoDict.subDict("CJ").lookupOrDefault("Rho", 0.0);
    CJ_MolWeight = thermoDict.subDict("CJ").lookupOrDefault("MolWeight", 0.0);
    CJ_Kminus = thermoDict.subDict("CJ").lookupOrDefault("Kminus", 0.0);
}

template<class ThermoType>
RDEMixtureBasic<ThermoType>::~RDEMixtureBasic()
{
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::CalculateCJ()
{
    bool isMolWeight;
    // Вычисления
    if(CJ_MolWeight > 1.0) // Если есть молярная масса
    {
        CJ_Rho = CJ_P*CJ_MolWeight/CJ_T/Foam::constant::thermodynamic::RR;
        isMolWeight = true;
    }else
    {
        CJ_MolWeight = CJ_Rho*CJ_T*Foam::constant::thermodynamic::RR/CJ_P;
        isMolWeight = false;
    }
    // Вычисляем Kminus, если его значение нулевое
    if(CJ_Kminus > Foam::SMALL);
    else
    {
        Foam::scalar up = CJ_Rho*Foam::exp(Edissociation/Foam::constant::thermodynamic::RR/CJ_T)*(1.0 - CJ_MolWeight/MolWeightMax)*(1.0 - CJ_MolWeight/MolWeightMax);
        Foam::scalar down = CJ_MolWeight*(CJ_MolWeight/MolWeightMin - 1.0)*Foam::pow(CJ_T/T0, beta/2.0)*Foam::pow((1.0 - Foam::exp(-theta/CJ_T)), beta);
        CJ_Kminus = up/down;
    }
    if(isMolWeight) Foam::Info << "CJ: (D = " << CJ_D << ", P = " << CJ_P/101325.0 << "atm, T = " << CJ_T << "K, MolWeight = " << CJ_MolWeight << ") Calculated: Rho = " << CJ_Rho << ", theta = " << theta << ", Kminus = " << CJ_Kminus << Foam::endl;
    else Foam::Info << "CJ: (D = " << CJ_D << ", P = " << CJ_P/101325.0 << "atm, T = " << CJ_T << "K, Rho = " << CJ_Rho << ") Calculated: MolWeight = " << CJ_MolWeight << ", theta = " << theta << ", Kminus = " << CJ_Kminus << Foam::endl;
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetMoles(const Foam::word & name) const
{
    if(!this->species().found(name)) return 0.0;
    if(!this->speciesData()(this->species()[name])) return 0.0;
    return this->speciesData()[this->species()[name]].moles();
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::read(const Foam::dictionary & thermoDict)
{
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::Inert(const Foam::label celli)
{
    Foam::scalar & Induction = Induction_[celli];
    Foam::scalar & FInduction = FInduction_[celli];
    Foam::scalar & MolWeight = MolWeight_[celli];
    Foam::scalar & FMolWeight = FMolWeight_[celli];

    FInduction = 0.0;
    FMolWeight = 0.0;
}

template<class ThermoType>
const word Foam::RDEMixtureBasic<ThermoType>::InfoElement(const Foam::label celli)
{
    Foam::scalar MolWeight = MolWeight_[celli];
    Foam::scalar Induction = Induction_[celli];
    return "Induction = " + Foam::name(Induction) + ", MolWeight = " + Foam::name(MolWeight);
}

template<class ThermoType>
const word Foam::RDEMixtureBasic<ThermoType>::InfoElement(const label patchi, const label facei)
{
    Foam::scalar MolWeight = MolWeight_.boundaryField()[patchi][facei];
    Foam::scalar Induction = Induction_.boundaryField()[patchi][facei];
    return "Induction = " + Foam::name(Induction) + ", MolWeight = " + Foam::name(MolWeight);
}

//--------------------------------------------------------------------------------------------------
template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetA(const Foam::scalar MolWeight, const Foam::scalar T) const
{
    if(T <= Tmin || T < 10.0) return 3.0/4.0*(MolWeight/MolWeightA + 1.0); // Если слишком малое T, => огромное theta/T
    if(theta/T < 0.001) return 3.0/4.0*(MolWeight/MolWeightA + 1.0) + 3.0/2.0*(MolWeight/MolWeightA - 1.0);
    return 3.0/4.0*(MolWeight/MolWeightA + 1.0) + 3.0/2.0*(MolWeight/MolWeightA - 1.0)*(theta/T)/(Foam::exp(theta/T) - 1.0);            // Николаев
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::Cv(const Foam::scalar P, const Foam::scalar T, const Foam::scalar MolWeight) const
{
    const Foam::scalar R = Foam::constant::thermodynamic::RR / MolWeight;
    return this->GetA(MolWeight, T) * R;
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetTmax(const Foam::scalar MolWeight, const Foam::scalar Uthermodinamical) const
{
    const Foam::scalar R = Foam::constant::thermodynamic::RR / MolWeight;
    const Foam::scalar A = 3.0/4.0*(MolWeight/this->MolWeightA + 1.0); //GetA при T->∞
    Foam::scalar Tmax = Uthermodinamical / A / R;
    Foam::scalar Tbound = Tinert;
    if(Tmax < Tbound)
    {
        Foam::Pout << "Tmax < " << Tbound << ", Told = " << Tmax;
        Tmax = Tbound;
        Foam::Pout << ", Tnew = " << Tmax << endl;
    }
    return Tmax;
}
#endif
