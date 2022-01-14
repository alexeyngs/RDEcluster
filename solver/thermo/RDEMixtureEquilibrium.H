#ifndef RDEMixtureEquilibrium_H
#define RDEMixtureEquilibrium_H

namespace Foam{
template<class ThermoType>
class RDEMixtureEquilibrium : public RDEMixtureBasic<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureEquilibrium(const RDEMixtureEquilibrium &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureEquilibrium &) = delete;
public:
    typedef ThermoType thermoType;
    void SetMixtureEquilibrium(const Foam::scalar Rho, const Foam::scalar E, const Foam::scalar Induction, const Foam::scalar alpha, const ThermoType & mixture)const;
    scalar MoleSpecieEquilibrium(const word name, const scalar Rho, const scalar E) const; // Равновесная термодинамика от (Rho, E)
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureEquilibrium(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureEquilibrium(){}
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureEquilibrium<" + ThermoType::typeName() + '>';
    }
    const ThermoType & Mixture(const ChemistryElement & Element) const;
    const ThermoType & cellMixture(const label celli) const;
    const ThermoType & cellT(scalar & Tout, const label celli, const scalar Hpotential, const scalar rho, const scalar Pav)const;
    const ThermoType & patchFaceMixture(const label patchi, const label facei) const;
    const ThermoType & patchFaceT(scalar & Tout, const label patchi, const label facei, const scalar Hpotential, const scalar rho, const scalar Pav) const;
    ChemistryElement GetFromHchemical(const ChemistryElement & old, scalar Hchemical)const;

    void SetWeightAndT(const ThermoType & mixture, const scalar Hpotential, const scalar rho, const scalar Pav, scalar & Induction, scalar & MolWeight, scalar & T)const;
    void Chemistry(const label celli, const scalar rho, const scalar T);
    void CorrectChemistry(ChemistryElement & Element);
    scalar InletButtMolWeight(const label patchi, const label facei, const ThermoType & mixture) const;
    scalar InletSideWallMolWeight(const label patchi, const label facei, const ThermoType & mixture) const;
    //--------------------------------------------------------------------------------------------------
    virtual scalar GetTmax(const scalar MolWeight, const scalar Uthermodinamical) const override;
};
} // End namespace Foam

template<class ThermoType>
Foam::RDEMixtureEquilibrium<ThermoType>::RDEMixtureEquilibrium(const dictionary& thermoDict, const fvMesh& mesh, const word& phaseName) : 
Foam::RDEMixtureBasic<ThermoType>(thermoDict, mesh, phaseName)
{
    const Foam::Time & runTime = mesh.time();
    this->MolWeightMax = GetMolWeightMaxVasilev<ThermoType>(this->species(), this->speciesData(), this->Sigma);
    if(this->Qmax == 0.0) this->Qmax = GetQMaxNikolaev<ThermoType>(this->species(), this->speciesData());
    // Энергия диссоциации от min до max:
    this->Edissociation = this->Qmax / this->Sigma / (1.0/this->MolWeightMin - 1.0/this->MolWeightMax);
    //-----------------------------------------------------------------
    Foam::Info << "Needs atoms(Oxigen) to total dissociate = " << this->fuel << Foam::endl;
    Foam::Info << "Qmax = " << this->Qmax << " J (" << this->Qmax / KiloCalorie << " kcal)" << Foam::endl;
    Foam::Info << "Edissociation = " << this->Edissociation/1000.0 << " J/mol (" << this->Edissociation/KiloCalorie/1000.0 << " kcal/mol)" <<  Foam::endl;
    Foam::Info << "Sigma = " << this->Sigma << Foam::endl;
    Foam::Info << "MolWeight0 = " << this->MolWeight0 << ", MolWeightA = " << this->MolWeightA << ", MolWeightMin = " << this->MolWeightMin << ", MolWeightMax = " << this->MolWeightMax << Foam::endl;
	Foam::Info << "CompareWeight = " << this->CompareWeight << Foam::endl;
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
            this->MolWeight_[i] = this->MolWeight0;
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
    this->Edissociation *= 1000.0;
    //----------------------------------------------------------------
    // Индукция
    Foam::scalar SumMoles = 0.0;
    this->Kinduction = 1.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > 0.0) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > 0.0)
    {
        scalar beta = this->speciesData()[i].moles() / SumMoles;
        this->Kinduction *= Foam::pow(this->speciesData()[i].Kinduction_, beta);
    }
    SumMoles = 0.0;
    this->Einduction = 0.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > 0.0) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > 0.0)
    {
        scalar beta = this->speciesData()[i].moles() / SumMoles;
        this->Einduction += this->speciesData()[i].Einduction_ * beta;
    }
    Foam::Info << "Einduction = " << this->Einduction << " J/mol, Kinduction = " << this->Kinduction << " mol*c/m3" << Foam::endl;
    this->Einduction *= 1000.0;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::Mixture(const ChemistryElement & Element) const
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
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::cellMixture(const label celli) const
{
    const Foam::IMixture * Mixture = dynamic_cast<const Foam::IMixture*>(this);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::cellMixture(celli);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = const_cast<Foam::volScalarField::FieldType &> (dynamic_cast<const Foam::volScalarField::FieldType &> (this->MolWeight_.primitiveField())) [celli];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_[celli];
    mixture.Mixture = const_cast<Foam::IMixture*>(Mixture);
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;

    const Foam::scalar RHO = this->Rho()[celli];
    const Foam::scalar HE = this->he()[celli];
    const Foam::scalar Sigma = mixture.W()*mixture.moles();           // Суммарная масса всей смеси
    const Foam::scalar MolWeightStar = GetMolWeightFromInductionVasilev(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
    Foam::scalar alpha = (MolWeight - MolWeightStar)/(this->CJ_MolWeight - MolWeightStar);
    if(alpha < 0.0) alpha = 0.0;
    if(alpha > 1.0) alpha = 1.0;
    SetMixtureEquilibrium(RHO, HE, Induction, alpha, mixture);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::cellT(scalar & Tout, const label celli, const scalar Hpotential, const scalar rho, const scalar Pav)const
{
    const ThermoType & mixture = cellMixture(celli);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = const_cast<Foam::volScalarField::FieldType &> (dynamic_cast<const Foam::volScalarField::FieldType &> (this->MolWeight_.primitiveField())) [celli];
    // Induction - индукция масса текущей ячейки
	Foam::scalar & Induction = const_cast<Foam::volScalarField::FieldType &> (dynamic_cast<const Foam::volScalarField::FieldType &> (this->Induction_.primitiveField())) [celli];
    //
    SetWeightAndT(mixture, Hpotential, rho, Pav, Induction, MolWeight, Tout);
    return mixture;
}


template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::patchFaceMixture(const Foam::label patchi, const Foam::label facei) const
{
    const Foam::IMixture * Mixture = dynamic_cast<const Foam::IMixture*>(this);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::patchFaceMixture(patchi, facei);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = const_cast<Foam::volScalarField::Boundary &> (dynamic_cast<const Foam::volScalarField::Boundary &> (this->MolWeight_.boundaryField())) [patchi][facei];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    mixture.Mixture = const_cast<Foam::IMixture*>(Mixture);
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
 
    const Foam::scalar RHO = this->Rho().boundaryField()[patchi][facei];
    const Foam::scalar HE = this->he().boundaryField()[patchi][facei];
    const Foam::scalar Sigma = mixture.W()*mixture.moles();           // Суммарная масса всей смеси
    const Foam::scalar MolWeightStar = GetMolWeightFromInductionVasilev(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
    Foam::scalar alpha = (MolWeight - MolWeightStar)/(this->CJ_MolWeight - MolWeightStar);
    if(alpha < 0.0) alpha = 0.0;
    if(alpha > 1.0) alpha = 1.0;
    SetMixtureEquilibrium(RHO, HE, Induction, alpha, mixture);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::patchFaceT(scalar & Tout, const Foam::label patchi, const Foam::label facei, const scalar Hpotential, const scalar rho, const scalar Pav) const
{
    const ThermoType & mixture = patchFaceMixture(patchi, facei);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = const_cast<Foam::volScalarField::Boundary &> (dynamic_cast<const Foam::volScalarField::Boundary &> (this->MolWeight_.boundaryField())) [patchi][facei];
    // Induction - индукция масса текущей ячейки
	Foam::scalar & Induction = const_cast<Foam::volScalarField::Boundary &> (dynamic_cast<const Foam::volScalarField::Boundary &> (this->Induction_.boundaryField())) [patchi][facei];
    //
    SetWeightAndT(mixture, Hpotential, rho, Pav, Induction, MolWeight, Tout);
    return mixture;
}

template<class ThermoType>
Foam::ChemistryElement Foam::RDEMixtureEquilibrium<ThermoType>::GetFromHchemical(const Foam::ChemistryElement & old, Foam::scalar Hchemical)const
{
    Foam::scalar Induction = std::get<0>(old);
    Foam::scalar MolWeight = std::get<1>(old);
    MolWeight = 1.0/(1.0/this->MolWeightMin + Hchemical/this->Edissociation);
    return std::make_tuple(Induction, MolWeight);
}

template<class ThermoType>
void Foam::RDEMixtureEquilibrium<ThermoType>::SetMixtureEquilibrium(const Foam::scalar Rho, const Foam::scalar E, const Foam::scalar Induction, const Foam::scalar alpha, const ThermoType & mixture)const
{
    // Суммарная масса смеси
    const Foam::scalar Sigma = mixture.W()*mixture.moles();
    if(Induction > SMALL)
    {
        mixture.MolWeight_ = GetMolWeightFromInductionVasilev(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
        mixture.Hchemical_ = GetHchemicalFromInductionVasilev(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
    }else
    {
        const Foam::scalar MolWeightStar = GetMolWeightFromInductionVasilev(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        const Foam::scalar HchemicalStar = GetHchemicalFromInductionVasilev(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        mixture.MolWeight_ = 0.0;       // молярная масса смеси ячейки
        mixture.Edissociation_ = 0.0;   // энергия диссоцияции ячейки
        Foam::scalar SummaMoles = 0.0;  // суммарное количество молей ячейки
        forAll(this->species(), i)
        {
            Foam::scalar mole = (1.0-alpha)*this->speciesData()[i].moles() + alpha*MoleSpecieEquilibrium(this->species()[i], Rho, E);
            SummaMoles += mole;  // Сумма всех молей
            mixture.MolWeight_ += this->speciesData()[i].MolWeight_ * mole;
            mixture.Edissociation_ += this->speciesData()[i].Edissociation_ * mole;
        }
        mixture.MolWeight_ /= SummaMoles;
        mixture.Edissociation_ /= SummaMoles;
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Foam::scalar Hdelta = HchemicalStar - mixture.Edissociation_*(1.0/MolWeightStar - 1.0/this->MolWeightMin);
        Foam::scalar Psi = (mixture.MolWeight_ < MolWeightStar) ? 1.0 : (mixture.MolWeight_ > this->CJ_MolWeight) ? 0.0 : (1.0/this->CJ_MolWeight - 1.0/mixture.MolWeight_) / (1.0/this->CJ_MolWeight - 1.0/MolWeightStar);
        mixture.Hchemical_ = mixture.Edissociation_*(1.0/mixture.MolWeight_ - 1.0/this->MolWeightMin) + Hdelta*Psi;
    }
}

template<class ThermoType>
inline void Foam::RDEMixtureEquilibrium<ThermoType>::SetWeightAndT(const ThermoType & mixture, const scalar Hpotential, const scalar rho, const scalar Pav, scalar & Induction, scalar & MolWeight, scalar & T)const
{
    Foam::scalar Hthermodinamical, Hchemical;
    mixture.MolWeight_ = MolWeight;
    Foam::scalar Tav = Pav*MolWeight/rho/Foam::constant::thermodynamic::RR;
    if(Tav < Tmin) Tav = Tmin;
    // Вычисление новой молярной массы
    Hthermodinamical = mixture.Hs(Pav, Tav); // Примерно определяем термодинамическую часть энтальпии
    //---------------------------------------------------------
    // Зная сумму и термодинамичекую часть энтальпии, определяемя химическую часть энтальпии
    Hchemical = Hpotential - Hthermodinamical;
    MolWeight = 1.0/(Hchemical/this->Edissociation + 1.0/this->MolWeightMin);
    if(MolWeight < this->MolWeightMin) MolWeight = this->MolWeightMin;
    if(MolWeight > this->MolWeightMax) MolWeight = this->MolWeightMax;
    mixture.MolWeight_ = MolWeight;
    // далее ищем температуру методом деления отрезка пополам
    Foam::scalar a,b;
    a = Tmin;
    b = this->GetTmax(Hthermodinamical, MolWeight);
    while(abs(a-b)/(a+b) > Foam::SMALL)
    {
        T = (a+b)/2.0;
        if(mixture.Hs(Pav, T) > Hthermodinamical) b = T;
        else a = T;
    }
    T = (a+b)/2.0;
    //---------------------------------------------------------
    Foam::Info << ", Hpotential = " << Hpotential << ", Hthermodinamical = " << Hthermodinamical << ", Hchemical = " << Hchemical;
}

template<class ThermoType>
void Foam::RDEMixtureEquilibrium<ThermoType>::Chemistry(const label celli, const scalar rho, const scalar T)
{
    Foam::scalar & Induction = this->Induction_[celli];
    Foam::scalar & FInduction = this->FInduction_[celli];
    Foam::scalar & MolWeight = this->MolWeight_[celli];
    Foam::scalar & FMolWeight = this->FMolWeight_[celli];
    const ThermoType & mixture = Foam::RDEMixtureBasic<ThermoType>::cellMixture(celli);
    const Foam::scalar fuel = mixture.needO()*mixture.moles();                         // Количество атомов кислорода для полного сгорания
    const Foam::scalar Sigma = mixture.W()*mixture.moles();                            // Суммарная масса всей смеси
    const Foam::scalar atoms = mixture.atoms()*mixture.moles();                        // Общее количество атомов с учетом молей
    //------------------------------------------------------
    // Стадия №1 - задержка воспламенения - уменьшение Induction от 1 до 0
    if(T > this->StartInduction)
    {
        Foam::scalar Concentration = 1.0;
        forAll(this->species(), i)
        {
            if(this->speciesData()[i].InductionPower() < Foam::SMALL) continue;
            Foam::scalar Ci = rho * this->speciesData()[i].moles() / Sigma * 1000.0;     // Концентрация компоненты смеси [кг/м3 * моль / кг] в моль/м3
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
            4.0*Kplus*this->CJ_Kminus*rho*(MolWeight/this->MolWeightMin - 1.0)*Foam::pow(T/T0, this->beta/2.0)*Foam::pow(1.0 - Foam::exp(-this->theta/T), this->beta)*Foam::exp(-this->Edissociation/Foam::constant::thermodynamic::RR/T);
        if(FWeight < this->CompareWeight) FWeight = this->CompareWeight;
        FMolWeight = rho*FWeight;
    }
}

template<class ThermoType>
void Foam::RDEMixtureEquilibrium<ThermoType>::CorrectChemistry(ChemistryElement & Element)
{
    Foam::scalar Induction = std::get<0>(Element);
    Foam::scalar MolWeight = std::get<1>(Element);
    //-------------------------------------
    if(Induction < 0.0) Induction = 0.0;
    if(Induction > Foam::SMALL)
    {
        MolWeight = this->MolWeight0;
    }else
    {    
        if(MolWeight < this->MolWeightMin) MolWeight = this->MolWeightMin;
        if(MolWeight > this->MolWeightMax) MolWeight = this->MolWeightMax;
    }
    //-------------------------------------
    std::get<0>(Element) = Induction;
    std::get<1>(Element) = MolWeight;
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureEquilibrium<ThermoType>::InletButtMolWeight(const label patchi, const label facei, const ThermoType & mixture) const
{
    return this->MolWeight0;
}


template<class ThermoType>
Foam::scalar Foam::RDEMixtureEquilibrium<ThermoType>::InletSideWallMolWeight(const label patchi, const label facei, const ThermoType & mixture) const
{
    return this->MolWeight0;
}

//--------------------------------------------------------------------------------------------------
template<class ThermoType>
Foam::scalar Foam::RDEMixtureEquilibrium<ThermoType>::GetTmax(const Foam::scalar MolWeight, const Foam::scalar Uthermodinamical) const
{
    const Foam::scalar R = Foam::constant::thermodynamic::RR / MolWeight;
    const Foam::scalar A = 3.0/4.0*(MolWeight/this->MolWeightA + 1.0); //GetA при T->∞
    Foam::scalar Tmax = Uthermodinamical / A / R;
    Foam::scalar Tbound = Tinert;
    if(Tmax < Tbound)
    {
        Info << "Tmax < " << Tbound << ", Told = " << Tmax;
        Tmax = Tbound;
        Info << ", Tnew = " << Tmax << endl;
    }
    return Tmax;
}

template<typename ThermoType>
Foam::scalar Foam::RDEMixtureEquilibrium<ThermoType>::MoleSpecieEquilibrium(const Foam::word name, const Foam::scalar Rho, const Foam::scalar E) const
{
    return 0.1;
}

#endif