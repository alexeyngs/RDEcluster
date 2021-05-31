#ifndef RDEMixtureExtension_H
#define RDEMixtureExtension_H
namespace Foam{
template<class ThermoType>
class RDEMixtureExtension : public RDEMixtureBasic<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureExtension(const RDEMixtureExtension &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureExtension &) = delete;
public:
    typedef ThermoType thermoType;
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureExtension(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureExtension();
    //- Return the instantiated type name
    static Foam::word typeName()
    {
        return "RDEMixtureExtension<" + ThermoType::typeName() + '>';
    }
    const ThermoType & Mixture(const ChemistryElement & Element) const;
    const ThermoType & cellMixture(const label celli) const;
    const ThermoType & cellT(scalar & Tout, const label celli, const scalar Hpotential, const scalar rho, const scalar Pav)const;
    const ThermoType & patchFaceMixture(const label patchi, const label facei) const;
    const ThermoType & patchFaceT(scalar & Tout, const label patchi, const label facei, const scalar Hpotential, const scalar rho, const scalar Pav) const;
    ChemistryElement GetFromHchemical(const ChemistryElement & old, scalar Hchemical)const;
    ChemistryElement CorrectElementAndT(const ThermoType & mixture, const ChemistryElement & old, const scalar Hpotential, const scalar rho, const scalar Pav, const scalar Tav, scalar & Tout)const;

    void SetWeightAndT(const ThermoType & mixture, const scalar Hpotential, const scalar rho, const scalar Pav, scalar molWeight, const scalar Induction, scalar & molWeightout, scalar & Tout)const;
    void Chemistry(const label celli, const scalar rho, const scalar T);
    void CorrectChemistry(ChemistryElement & Element)const;
    scalar InletButtMolWeight(const label patchi, const label facei, const ThermoType & mixture) const;
    scalar InletSideWallMolWeight(const label patchi, const label facei, const ThermoType & mixture) const;
    //--------------------------------------------------------------------------------------------------
    virtual scalar GetTmax(const scalar MolWeight, const scalar Uthermodinamical) const override;
    // Дополнительные специфические поля и методы
    scalar Hchemical(const scalar Induction, scalar MolWeight)const;
    mutable scalar MolWeightStar;
    mutable scalar HchemicalStar;
    mutable scalar Hdelta;
    //--------------------------------------------------------------------------------------------------
    virtual scalar GetMolWeightFromInduction(const Foam::scalar Induction) const;
    virtual scalar GetHchemicalFromInduction(const Foam::scalar Induction) const;
    virtual scalar GetMolWeightMax(const Foam::dictionary & thermoDictionary) const;
    virtual scalar GetQMax(const Foam::dictionary & thermoDictionary) const;
    virtual scalar GetTheta() const;
};
} // End namespace Foam

template<class ThermoType>
Foam::RDEMixtureExtension<ThermoType>::RDEMixtureExtension(const Foam::dictionary & thermoDictionary, const Foam::fvMesh & mesh, const Foam::word & phaseName):
Foam::RDEMixtureBasic<ThermoType>(thermoDictionary, mesh, phaseName)
{
    const Foam::Time & runTime = mesh.time();
    // На данный момент известны Sigma и компоненты смеси
    this->MolWeight0 = GetMolWeightFromInduction(1.0);
    this->MolWeightStar = GetMolWeightFromInduction(0.0);
    this->MolWeightMax = GetMolWeightMax(thermoDictionary);
    if(this->Qmax == 0.0) this->Qmax = GetQMax(thermoDictionary);
    // Энергия диссоциации от * до max:
//    this->Edissociation = (Qmax - Qstar) / Sigma / (1.0/MolWeightStar - 1.0/MolWeightMax);
    // Энергия диссоциации от min до max:
    this->Edissociation = this->Qmax / this->Sigma / (1.0/this->MolWeightMin - 1.0/this->MolWeightMax);
    this->HchemicalStar = GetHchemicalFromInduction(0.0);
    this->Edissociation *= 1000.0; //С этого момента берем Дж/киломоль
    this->Hdelta = this->HchemicalStar - this->Edissociation*(1.0/this->MolWeightStar - 1.0/this->MolWeightMin);
    if(this->theta < 1.0) this->theta = GetTheta();
    //-----------------------------------------------------------------
    Foam::Info << "Needs atoms(Oxigen) to total dissociate = " << this->fuel << Foam::endl;
    Foam::Info << "Qmax = " << this->Qmax << " J (" << this->Qmax / KiloCalorie << " kcal)" << Foam::endl;
    Foam::Info << "Edissociation = " << this->Edissociation/1000.0 << " J/mol (" << this->Edissociation/KiloCalorie/1000.0 << " kcal/mol)" <<  Foam::endl;
    Foam::Info << "Sigma = " << this->Sigma << ", Hchemical* = "<< this->HchemicalStar << ", Hdelta = " << this->Hdelta << Foam::endl;
    Foam::Info << "MolWeight0 = " << this->MolWeight0 << ", MolWeight* = " << this->MolWeightStar << ", MolWeightA = " << this->MolWeightA << ", MolWeightMin = " << this->MolWeightMin << ", MolWeightMax = " << this->MolWeightMax << Foam::endl;
    Foam::RDEMixtureBasic<ThermoType>::CalculateCJ();
    //-----------------------------------------------------------------
    forAll(multiComponentMixture<ThermoType>::speciesData(), i)
    {
        multiComponentMixture<ThermoType>::setActive(i);
    }

    // Заполнение поля молярной массы, но только в том случае, когда расчет начинается с 0
    if(runTime.value() < Foam::SMALL)
    {
        // Заполнение поля в центре области
        forAll(basicMultiComponentMixture::Y(0) , i)
        {
            const ThermoType & mixture = multiComponentMixture<ThermoType>::cellMixture(i);
            Foam::scalar Induction = this->Induction_[i];
            this->MolWeight_[i] = GetMolWeightFromInduction(Induction);
        }
        // Заполнение поля на границах, кроме определенной границы "inlet"
        forAll(this->MolWeight_.boundaryField(), patchi)
        {
            if(this->MolWeight_.boundaryField()[patchi].patch().name() == "inlet") continue;
            forAll(this->MolWeight_.boundaryField()[patchi], facei)
            {
                this->MolWeight_.boundaryFieldRef()[patchi][facei] = this->MolWeight0;
            }
        }
    }

    // Гарантированное заполнение поля границы "inlet", вне зависимости от времени старта программы
    forAll(this->MolWeight_.boundaryField(), patchi)
    {
        if(this->MolWeight_.boundaryField()[patchi].patch().name() == "inlet")
            forAll(this->MolWeight_.boundaryField()[patchi], facei)
            {
                const ThermoType & mixture = patchFaceMixture(patchi, facei);
                this->MolWeight_.boundaryFieldRef()[patchi][facei] = InletButtMolWeight(patchi, facei, mixture);
            }
    }
    //----------------------------------------------------------------
    // Индукция
    Foam::scalar SumMoles = 0.0;
    this->Kinduction = 1.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > Foam::SMALL) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > Foam::SMALL)
    {
        Foam::scalar beta = this->speciesData()[i].moles() / SumMoles;
        this->Kinduction *= Foam::pow(this->speciesData()[i].Kinduction_, beta);
    }
    SumMoles = 0.0;
    this->Einduction = 0.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > Foam::SMALL) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > Foam::SMALL)
    {
        Foam::scalar beta = this->speciesData()[i].moles() / SumMoles;
        this->Einduction += this->speciesData()[i].Einduction_ * beta;
    }
    Foam::Info << "Einduction = " << this->Einduction << " J/mol, Kinduction = " << this->Kinduction << " mol*c/m3" << Foam::endl;
    this->Einduction *= 1000.0;
}

template<class ThermoType>
Foam::RDEMixtureExtension<ThermoType>::~RDEMixtureExtension()
{
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureExtension<ThermoType>::Mixture(const ChemistryElement & Element) const
{
    Foam::scalar Induction = std::get<0>(Element);
    Foam::scalar MolWeight = std::get<1>(Element);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::cellMixture(0);
    mixture.Mixture = const_cast<Foam::IMixture*>(dynamic_cast<const Foam::IMixture*>(this));
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
    mixture.MolWeight_ = MolWeight;
    mixture.Hchemical_ = Hchemical(Induction, MolWeight);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureExtension<ThermoType>::cellMixture(const label celli) const
{
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::cellMixture(celli);
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = dynamic_cast<const Foam::volScalarField::FieldType &> (this->Induction_.primitiveField()) [celli];
    mixture.Mixture = const_cast<Foam::IMixture*>(dynamic_cast<const Foam::IMixture*>(this));
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
    mixture.MolWeight_ = this->MolWeight_.primitiveField() [celli];
    mixture.Hchemical_ = Hchemical(Induction, mixture.MolWeight_);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureExtension<ThermoType>::patchFaceMixture(const Foam::label patchi, const Foam::label facei) const
{
    const Foam::IMixture * Mixture = dynamic_cast<const Foam::IMixture*>(this);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::patchFaceMixture(patchi, facei);
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    mixture.Mixture = const_cast<Foam::IMixture*>(Mixture);
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
    mixture.MolWeight_ = this->MolWeight_.boundaryField() [patchi][facei];
    mixture.Hchemical_ = Hchemical(Induction, mixture.MolWeight_);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureExtension<ThermoType>::cellT(scalar & Tout, const label celli, const scalar Hpotential, const scalar rho, const scalar Pav)const
{
    const ThermoType & mixture = cellMixture(celli);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & molWeight = const_cast<Foam::volScalarField::FieldType &> (dynamic_cast<const Foam::volScalarField::FieldType &> (this->MolWeight_.primitiveField())) [celli];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_[celli];
    //
    Foam::scalar molWeightout;
    SetWeightAndT(mixture, Hpotential, rho, Pav, molWeight, Induction, molWeightout, Tout);
    molWeight = molWeightout;
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureExtension<ThermoType>::patchFaceT(scalar & Tout, const Foam::label patchi, const Foam::label facei, const scalar Hpotential, const scalar rho, const scalar Pav) const
{
    const ThermoType & mixture = patchFaceMixture(patchi, facei);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & molWeight = const_cast<Foam::volScalarField::Boundary &> (dynamic_cast<const Foam::volScalarField::Boundary &> (this->MolWeight_.boundaryField())) [patchi][facei];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    //
    Foam::scalar molWeightout;
    SetWeightAndT(mixture, Hpotential, rho, Pav, molWeight, Induction, molWeightout, Tout);
    molWeight = molWeightout;
    return mixture;
}


template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::Hchemical(const Foam::scalar Induction, Foam::scalar MolWeight)const
{
const Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
return this->Edissociation*(1.0/MolWeight - 1.0/this->MolWeightMin) + Hdelta*Psi;
//-------------------------------------------
    if(Induction > Foam::SMALL)
    {
	    const Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
        return GetHchemicalFromInduction(Induction);
    }else
    {
        const Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
        return this->Edissociation*(1.0/MolWeight - 1.0/this->MolWeightMin) + Hdelta*Psi;
    }
}

template<class ThermoType>
Foam::ChemistryElement Foam::RDEMixtureExtension<ThermoType>::GetFromHchemical(const Foam::ChemistryElement & old, Foam::scalar Hchemical)const
{
    Foam::scalar Induction = std::get<0>(old);
    Foam::scalar MolWeight = std::get<1>(old);
    // Если старая индукция положительная, то пробуем получить новую индукцию из химической части энтальпии
    /*if(Induction > Foam::SMALL)
    {
        Induction = GetInductionFromHchemical(Hchemical);
        if(Induction < 0.0) Induction = 0.0;
        if(Induction > 1.0) Induction = 1.0;
        // Теперь от 0 до 1
        if(Induction > Foam::SMALL) MolWeight = GetMolWeightFromInduction(Induction);
        else
        {
            Induction = 0.0;
            const Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
            // Вычисление новой молярной массы
            MolWeight = 1.0/((Hchemical - this->Hdelta*Psi)/this->Edissociation + 1.0/this->MolWeightMin);
        }
        
    }else*/
    {
        Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
        // Вычисление новой молярной массы
        MolWeight = 1.0/((Hchemical - this->Hdelta*Psi)/this->Edissociation + 1.0/this->MolWeightMin);
        // дубль...
        Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/MolWeight) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
        MolWeight = 1.0/((Hchemical - this->Hdelta*Psi)/this->Edissociation + 1.0/this->MolWeightMin);
    }
    //if(MolWeight < this->MolWeightMin || MolWeight > this->MolWeightMax) return old;
    return std::make_tuple(Induction, MolWeight);
}


template<class ThermoType>
inline void Foam::RDEMixtureExtension<ThermoType>::SetWeightAndT(const ThermoType & mixture, const scalar Hpotential, const scalar rho, const scalar Pav, scalar molWeight, const scalar Induction, scalar & molWeightout, scalar & Tout)const
{
    Foam::scalar Hthermodinamical, Hchemical;
    //if(molWeight < MolWeightStar) molWeight = MolWeightStar;
    mixture.MolWeight_ = molWeight;
    Foam::scalar Tav = Pav*molWeight/rho/Foam::constant::thermodynamic::RR;
    if(Tav < Tmin) Tav = Tmin;
    // Вычисление новой молярной массы
    /*if(Induction > Foam::SMALL)
    {
        molWeight = GetMolWeightFromInduction(Induction);
        Hchemical = GetHchemicalFromInduction(Induction);
        Hthermodinamical = Hpotential - Hchemical;
        if(Hthermodinamical < 0.0) Hthermodinamical = mixture.Hs(Pav, Tav); // Примерно определяем термодинамическую часть энтальпии
    }else*/
    {
        Hthermodinamical = mixture.Hs(Pav, Tav); // Примерно определяем термодинамическую часть энтальпии
        //---------------------------------------------------------
        // Зная сумму и термодинамичекую часть энтальпии, определяемя химическую часть энтальпии
        Hchemical = Hpotential - Hthermodinamical;
        ChemistryElement old = std::make_tuple(Induction, molWeight);
        ChemistryElement element = GetFromHchemical(old, Hchemical);
        CorrectChemistry(element);
        molWeight = std::get<1>(element);
        // Дубль...
        Tav = Pav*molWeight/rho/Foam::constant::thermodynamic::RR;
        Hthermodinamical = mixture.Hs(Pav, Tav); // Примерно определяем термодинамическую часть энтальпии
        Hchemical = Hpotential - Hthermodinamical;
        old = element;
        element = GetFromHchemical(old, Hchemical);
        CorrectChemistry(element);
        molWeight = std::get<1>(element);
    }
    mixture.MolWeight_ = molWeight;
    // далее ищем температуру методом деления отрезка пополам
    Foam::scalar a,b,T;
    a = Tmin;
    b = this->GetTmax(Hthermodinamical, molWeight);
    while(abs(a-b)/(a+b) > Foam::SMALL)
    {
        T = (a+b)/2.0;
        if(mixture.Hs(Pav, T) > Hthermodinamical) b = T;
        else a = T;
    }
    T = (a+b)/2.0;
    //---------------------------------------------------------
    Foam::Info << ", Hpotential = " << Hpotential << ", Hthermodinamical = " << Hthermodinamical << ", Hchemical = " << Hchemical;
    molWeightout = molWeight;
    Tout = T;
}


template<class ThermoType>
ChemistryElement Foam::RDEMixtureExtension<ThermoType>::CorrectElementAndT(const ThermoType & mixture, const ChemistryElement & old, const scalar Hpotential, const scalar rho, const scalar Pav, const scalar Tav, scalar & Tout)const
{
    ChemistryElement Elementout = old;
    Foam::scalar Induction = std::get<0>(old);
    Foam::scalar MolWeight = std::get<1>(old);
    Foam::scalar Hthermodinamical, Hchemical;
    mixture.MolWeight_ = MolWeight;
    // Вычисление новой молярной массы
    /*if(Induction > Foam::SMALL)
    {
        MolWeight = GetMolWeightFromInduction(Induction);
        Hchemical = GetHchemicalFromInduction(Induction);
        Hthermodinamical = Hpotential - Hchemical;
        // Примерно определяем термодинамическую часть энтальпии
        if(Hthermodinamical < 0.0) Hthermodinamical = mixture.Hs(Pav, Tav);
    }else*/
    {
        // Примерно определяем термодинамическую часть энтальпии
        Hthermodinamical = mixture.Hs(Pav, Tav);
        // Зная сумму и термодинамичекую часть энтальпии, определяемя химическую часть энтальпии
        Hchemical = Hpotential - Hthermodinamical;
        Elementout = GetFromHchemical(old, Hchemical);
        CorrectChemistry(Elementout);
        MolWeight = std::get<1>(Elementout);
        // Дубль...
        Tav = Pav*MolWeight/rho/Foam::constant::thermodynamic::RR;
        // Примерно определяем термодинамическую часть энтальпии
        Hthermodinamical = mixture.Hs(Pav, Tav);
        Hchemical = Hpotential - Hthermodinamical;
        Elementout = GetFromHchemical(Elementout, Hchemical);
        CorrectChemistry(Elementout);
        MolWeight = std::get<1>(Elementout);
    }
    mixture.MolWeight_ = MolWeight;
    // далее ищем температуру методом деления отрезка пополам
    Foam::scalar a,b,T;
    a = Tmin;
    b = this->GetTmax(Hthermodinamical, MolWeight);
    while(abs(a-b)/(a+b) > Foam::SMALL)
    {
        T = (a+b)/2.0;
        if(mixture.Hs(Pav, T) > Hthermodinamical) b = T;
        else a = T;
    }
    T = (a+b)/2.0;
    Tout = T;
    std::get<1>(Elementout) = MolWeight;
    return Elementout;
}

template<class ThermoType>
void Foam::RDEMixtureExtension<ThermoType>::Chemistry(const label celli, const scalar rho, const scalar T)
{
    Foam::scalar & Induction = this->Induction_[celli];
    Foam::scalar & FInduction = this->FInduction_[celli];
    Foam::scalar & MolWeight = this->MolWeight_[celli];
    Foam::scalar & FMolWeight = this->FMolWeight_[celli];
    //------------------------------------------------------
    // Стадия №1 - задержка воспламенения - уменьшение Induction от 1 до 0
    if(T > this->StartInduction)
    {
        Foam::scalar Concentration = 1.0;
        forAll(this->species(), i)
        {
            if(this->speciesData()[i].InductionPower() < Foam::SMALL) continue;
            Foam::scalar Ci = rho * this->speciesData()[i].moles() / this->Sigma * 1000.0;     // Концентрация компоненты смеси [кг/м3 * моль / кг] в моль/м3
            Concentration *= Foam::pow(Ci, this->speciesData()[i].InductionPower());
        }
        Foam::scalar tind = this->Kinduction / Concentration * Foam::exp(this->Einduction / Foam::constant::thermodynamic::RR / T);
        FInduction = -1.0/tind;
    }else
    {
        FInduction = 0.0;
    }

    // Стадия №2 - химическая реакция
    if(Induction > Foam::SMALL || T < Tinert)
    {
        FMolWeight = 0.0;
    }
    else    // Вычисление приращения молекулярной массы
    {
        Foam::scalar FWeight =
            4.0*Kplus*rho*rho/MolWeight * (1.0 - MolWeight/this->MolWeightMax)*(1.0 - MolWeight/this->MolWeightMax) -
            4.0*Kplus*this->CJ_Kminus*rho*(MolWeight/this->MolWeightMin - 1.0)*Foam::pow(T/T0, this->beta/2.0)*Foam::pow(1.0 - Foam::exp(-this->theta/T), Foam::RDEMixtureBasic<ThermoType>::beta)*Foam::exp(-this->Edissociation/Foam::constant::thermodynamic::RR/T);        
        if(FWeight < -10.0) FWeight = -10.0;
        FMolWeight = rho*FWeight;
    }
}

template<class ThermoType>
void Foam::RDEMixtureExtension<ThermoType>::CorrectChemistry(ChemistryElement & Element)const
{
    Foam::scalar Induction = std::get<0>(Element);
    Foam::scalar MolWeight = std::get<1>(Element);
    //-------------------------------------
    if(MolWeight < this->MolWeightMin) MolWeight = this->MolWeightMin;
    if(MolWeight > this->MolWeightMax) MolWeight = this->MolWeightMax;
//    if(Induction < 0.0) Induction = 0.0;
//    if(Induction > Foam::SMALL)
//    {
//        //MolWeight = GetMolWeightFromInduction(Induction);
//    }else
//    {
//        if(MolWeight < this->MolWeightMin) MolWeight = this->MolWeightMin;
//        if(MolWeight > this->MolWeightMax) MolWeight = this->MolWeightMax;
//    }
    //-------------------------------------
    std::get<0>(Element) = Induction;
    std::get<1>(Element) = MolWeight;
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::InletButtMolWeight(const label patchi, const label facei, const ThermoType & mixture) const
{
    Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    Foam::scalar MolWeight = GetMolWeightFromInduction(Induction);
    return MolWeight;
}


template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::InletSideWallMolWeight(const label patchi, const label facei, const ThermoType & mixture) const
{
    Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    Foam::scalar MolWeight = GetMolWeightFromInduction(Induction);
    return MolWeight;
}

//--------------------------------------------------------------------------------------------------
template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetTmax(const Foam::scalar MolWeight, const Foam::scalar Uthermodinamical) const
{
    const Foam::scalar R = Foam::constant::thermodynamic::RR / MolWeight;
    const Foam::scalar A = 3.0/4.0*(MolWeight/this->MolWeightA + 1.0); //GetA при T->∞
    Foam::scalar Tmax = Uthermodinamical / A / R;
    Foam::scalar Tbound = Tinert;
    if(Tmax < Tbound-10)
    {
        Info << "Tmax < " << Tbound << ", Told = " << Tmax << ", MolWeight = " << MolWeight;
        Tmax = Tbound;
        Info << ", Tnew = " << Tmax << endl;
    }
    return Tmax;
}

//--------------------------------------------------------------------------------------------------
template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetMolWeightFromInduction(const Foam::scalar Induction) const
{
    return GetMolWeightFromInductionVasilev<ThermoType>(this->species(), this->speciesData(), Induction, this->powerN, this->Sigma);
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetHchemicalFromInduction(const Foam::scalar Induction) const
{
    return GetHchemicalFromInductionVasilev(this->species(), this->speciesData(), Induction, this->powerN, this->Sigma);
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetMolWeightMax(const Foam::dictionary & thermoDictionary) const
{
    return GetMolWeightMaxVasilev<ThermoType>(this->species(), this->speciesData(), this->Sigma);
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetQMax(const Foam::dictionary & thermoDictionary) const
{
    return GetQMaxVasilev<ThermoType>(this->species(), this->speciesData());
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureExtension<ThermoType>::GetTheta() const
{
    return GetThetaVasilev(this->species(), this->speciesData());
}

#endif