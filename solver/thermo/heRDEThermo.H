#ifndef heRDEThermo_H
#define heRDEThermo_H

#include "rhoRDEThermo.H"
#include "RDEThermo.H"
#include "turbulentFluidThermoModel.H"

extern Foam::scalar Tmin;
//const Foam::scalar RhoLimit = 6.0;

namespace Foam
{
template<class BasicPsiThermo, class MixtureType>
class heRDEThermo : public heThermo<BasicPsiThermo, MixtureType>
{
    //- Calculate the thermo variables
    void calculateTHE();
    void calculategetT();
    //- Construct as copy (not implemented)
    heRDEThermo(const heRDEThermo<BasicPsiThermo, MixtureType>&);
public:
    //- Runtime type information
    TypeName("heRDEThermo");
    // Constructors
    //- Construct from mesh and phase name
    heRDEThermo(const fvMesh &, const word & phaseName);
    //- Destructor
    virtual ~heRDEThermo();
    // Member functions
    virtual void correct();
    // от класса rhoRDEThermo
    virtual volScalarField & CorrectRho();
    virtual const volScalarField & Rho() const override;
    virtual volScalarField & MolWeight() override;
    virtual volScalarField & FMolWeight() override;
    virtual volScalarField & Induction() override;
    virtual volScalarField & FInduction() override;
    virtual void chemistry(const Foam::scalar CameraLength) override;
    virtual void CorrectChemistry() override;
    virtual void CorrectChemistry(ChemistryElement & Element) override;
    virtual scalar Dcj() override;
    virtual void CorrectErrors(bool UseChemistry, volScalarField & rho) override;
    virtual ChemistryElement GetChemistry(label celli) override;
    virtual void SetChemistry(label celli, const ChemistryElement & Element) override;
    virtual ChemistryElement GetChemistryFromHchemical(const ChemistryElement & old, scalar Hchemical) override;
    virtual ChemistryElement GetChemistryFromHpotential(scalar rho, scalar P, scalar gamma, scalar Hpotential, const ChemistryElement & old) override;
    virtual scalar Hthermodynamical(label celli) override;
    virtual scalar Hchemical(label celli) override;
    virtual scalar Hthermodynamical(scalar P, scalar rho, const ChemistryElement & Element) override;
    virtual scalar Hchemical(const ChemistryElement & Element) override;
    virtual scalar Hpotential(scalar P, scalar pho, const ChemistryElement & Element) override;
    virtual scalar MolWeightStar() override;
    virtual scalar patchInletButtMolWeight(const label patchi, const label facei) override;
    virtual scalar patchInletSideWallMolWeight(const label patchi, const label facei) override;
    virtual scalar patchOutletMolWeight(const label patchi, const label facei) override;
    virtual const word InfoElement(const label celli) override;
    virtual const word InfoElement(const label patchi, const label facei) override;
};
} // End namespace Foam//

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::calculateTHE()
{
    const Foam::fvMesh & mesh = this->p_.mesh();
    Foam::scalarField & hCells = this->he();
    Foam::scalarField & pcell = this->p_.primitiveFieldRef();
    Foam::scalarField & Tcell = this->T_.primitiveFieldRef();
    Foam::scalarField & rhocell = this->rho_.primitiveFieldRef();
    Foam::scalarField & psiCells = this->psi_.primitiveFieldRef();    
    Foam::scalarField & muCells = this->mu_.primitiveFieldRef();
    Foam::scalarField & alphaCells = this->alpha_.primitiveFieldRef();
    forAll(Tcell, celli)
    {
        typename MixtureType::thermoType & mixture = const_cast<typename MixtureType::thermoType &>(this->cellMixture(celli));
        const Foam::point & Position = mesh.C()[celli];
        Foam::scalar Uchemical = mixture.Hc();
        Foam::scalar Uthermodynamical = hCells[celli] - Uchemical;
        if (Uthermodynamical > 0.0)
        {
            Tcell[celli] = mixture.THE(hCells[celli], pcell[celli], Tcell[celli]);
            if(Tcell[celli] < Tmin)
            {
                Foam::Pout << "(old THE)T < " << Tmin << "K!, p = " << pcell[celli]/101325.0 << "atm, T = " << Tcell[celli] << "K, rho = " << rhocell[celli] << ", mu = " << mixture.W() << ", Position = " << Position << ", Uthermodynamical = " << Uthermodynamical << ", Uchemical = " << Uchemical << endl;
                //<<<
                pcell[celli] = this->Paverage[celli];
                if(rhocell[celli] < 0.0) this->rho_[celli] = this->Rhoaverage[celli];
//if(rhocell[celli] > RhoLimit) rhocell[celli] = RhoLimit;
                Tcell[celli] = pcell[celli] / this->rho_[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
                Uthermodynamical = this->p_[celli] / this->rho_[celli] * mixture.Mixture->GetA(mixture.W(), Tcell[celli]);
                Uchemical = hCells[celli] - Uthermodynamical;
                // Получение нового химического элемента
                std::tie(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]), Uchemical);
                mixture = this->cellMixture(celli);
                //!!!
                Tcell[celli] = pcell[celli] / rhocell[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
                Uthermodynamical = pcell[celli] / rhocell[celli] * mixture.Mixture->GetA(mixture.W(), Tcell[celli]);
                Uchemical = mixture.Hc();
                hCells[celli] = Uthermodynamical + Uchemical;
                //>>>
                Tcell[celli] = mixture.THE(hCells[celli], pcell[celli], Tcell[celli]);
                Foam::Pout << "(new THE)T < " << Tmin << "K!, p = " << pcell[celli]/101325.0 << "atm, T = " << Tcell[celli] << "K, rho = " << rhocell[celli] << ", mu = " << mixture.W() << ", Position = " << Position << ", Uthermodynamical = " << Uthermodynamical << ", Uchemical = " << Uchemical << endl;
            }
        }else
        {
            Foam::Pout << "(old THE U<0)T < " << Tmin << "K!, p = " << pcell[celli]/101325.0 << "atm, T = " << Tcell[celli] << "K, rho = " << rhocell[celli] << ", Induction = " << MixtureType::Induction_[celli] << ", MolWeight = " << MixtureType::MolWeight_[celli] << ", Position = " << Position << ", Uthermodynamical = " << Uthermodynamical << ", Uchemical = " << Uchemical << endl;
            //<<<
            pcell[celli] = this->Paverage[celli];
            if(rhocell[celli] < 0.0) rhocell[celli] = this->Rhoaverage[celli];
//if(rhocell[celli] > RhoLimit) rhocell[celli] = RhoLimit;
            Tcell[celli] = pcell[celli] / rhocell[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
            Uthermodynamical = pcell[celli] / rhocell[celli] * mixture.Mixture->GetA(mixture.W(), Tcell[celli]);
            Uchemical = hCells[celli] - Uthermodynamical;
            // Получение нового химического элемента
            std::tie(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]), Uchemical);
            mixture = this->cellMixture(celli);
            //!!!
            Tcell[celli] = pcell[celli] / rhocell[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
            Uthermodynamical = pcell[celli] / rhocell[celli] * mixture.Mixture->GetA(mixture.W(), Tcell[celli]);
            Uchemical = mixture.Hc();
            hCells[celli] = Uthermodynamical + Uchemical;
            //>>>
            Tcell[celli] = mixture.THE(hCells[celli], pcell[celli], Tcell[celli]);
            Foam::Pout << "(new THE U<0)T < " << Tmin << "K!, p = " << pcell[celli]/101325.0 << "atm, T = " << Tcell[celli] << "K, rho = " << rhocell[celli] << ", Induction = " << MixtureType::Induction_[celli] << ", MolWeight = " << MixtureType::MolWeight_[celli] << ", Position = " << Position << ", Uthermodynamical = " << Uthermodynamical << ", Uchemical = " << Uchemical << endl;
        }
        // Заполнение остальных параметров
        psiCells[celli] = mixture.psi(pcell[celli], Tcell[celli]);
        rhocell[celli] = mixture.rho(pcell[celli], Tcell[celli]);
        muCells[celli] = mixture.mu(pcell[celli], Tcell[celli]);
        alphaCells[celli] = mixture.alphah(pcell[celli], Tcell[celli]);
    }

    // Граничные условия
    Foam::volScalarField::Boundary & pBoundaryField = this->p_.boundaryFieldRef();
    Foam::volScalarField::Boundary & TBoundaryField = this->T_.boundaryFieldRef();
    Foam::volScalarField::Boundary & psiBoundaryField = this->psi_.boundaryFieldRef();
    Foam::volScalarField::Boundary & rhoBoundaryField = this->rho_.boundaryFieldRef();
    Foam::volScalarField::Boundary & heBoundaryField = this->he().boundaryFieldRef();
    Foam::volScalarField::Boundary & muBoundaryField = this->mu_.boundaryFieldRef();
    Foam::volScalarField::Boundary & alphaBoundaryField = this->alpha_.boundaryFieldRef();
    forAll(this->T_.boundaryField(), patchi)
    {
        Foam::fvPatchScalarField & BoundaryP = pBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryT = TBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryPsi = psiBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryRho = rhoBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryHE = heBoundaryField[patchi];
        Foam::fvPatchScalarField & pmu = muBoundaryField[patchi];
        Foam::fvPatchScalarField & palpha = alphaBoundaryField[patchi];

        /*if (BoundaryT.fixesValue())
        {
            forAll(BoundaryT, facei)
            {
                const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
                BoundaryHE[facei] = mixture.HE(BoundaryP[facei], BoundaryT[facei]);
                ppsi[facei] = mixture.psi(BoundaryP[facei], BoundaryT[facei]);
                prho[facei] = mixture.rho(BoundaryP[facei], BoundaryT[facei]);
                pmu[facei] = mixture.mu(BoundaryP[facei], BoundaryT[facei]);
                palpha[facei] = mixture.alphah(BoundaryP[facei], BoundaryT[facei]);
            }
        }
        else*/
        {
            forAll(BoundaryT, facei)
            {
                typename MixtureType::thermoType & mixture = const_cast<typename MixtureType::thermoType &>(this->patchFaceMixture(patchi, facei));
                if(BoundaryHE[facei] - mixture.Hc() > 0.0) // Uthermodynamical
                {
                    BoundaryT[facei] = mixture.THE(BoundaryHE[facei], BoundaryP[facei], BoundaryT[facei]);
                    if(BoundaryT[facei] < Tmin)
                    {
                        const Foam::point & Position = mesh.boundary()[patchi].Cf()[facei];
                        Foam::Pout << "(old boundary THE)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, mu = " << mixture.W() << ", Position = " << Position << endl;
                        //<<<
                        BoundaryP[facei] = this->Paverage.boundaryField()[patchi][facei];
                        if(BoundaryRho[facei] < 0.0) BoundaryRho[facei] = this->Rhoaverage.boundaryField()[patchi][facei];
//if(BoundaryRho[facei] > RhoLimit) BoundaryRho[facei] = RhoLimit;
                        BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                        Foam::scalar Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                        Foam::scalar Uchemical = BoundaryHE[facei] - Uthermodynamical;
                        // Получение нового химического элемента
                        std::tie(MixtureType::Induction_.boundaryFieldRef()[patchi][facei], MixtureType::MolWeight_.boundaryFieldRef()[patchi][facei]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_.boundaryField()[patchi][facei], MixtureType::MolWeight_.boundaryField()[patchi][facei]), Uchemical);
                        mixture = this->patchFaceMixture(patchi, facei);
                        //!!!
                        BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                        Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                        Uchemical = mixture.Hc();
                        BoundaryHE[facei] = Uthermodynamical + Uchemical;
                        //>>>
                        BoundaryT[facei] = mixture.THE(BoundaryHE[facei], BoundaryP[facei], BoundaryT[facei]);
                        Foam::Pout << "(new boundary THE)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, mu = " << mixture.W() << ", Position = " << Position << endl;
                    }
                }else
                {
                    const Foam::point & Position = mesh.boundary()[patchi].Cf()[facei];
                    Foam::Pout << "(old boundary THE U<0)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, rho = " << BoundaryRho[facei] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
                    //<<<
                    BoundaryP[facei] = this->Paverage.boundaryField()[patchi][facei];
                    if(BoundaryRho[facei] < 0.0) BoundaryRho[facei] = this->Rhoaverage.boundaryField()[patchi][facei];
//if(BoundaryRho[facei] > RhoLimit) BoundaryRho[facei] = RhoLimit;
                    BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                    Foam::scalar Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                    Foam::scalar Uchemical = BoundaryHE[facei] - Uthermodynamical;
                    // Получение нового химического элемента
                    std::tie(MixtureType::Induction_.boundaryFieldRef()[patchi][facei], MixtureType::MolWeight_.boundaryFieldRef()[patchi][facei]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_.boundaryField()[patchi][facei], MixtureType::MolWeight_.boundaryField()[patchi][facei]), Uchemical);
                    mixture = this->patchFaceMixture(patchi, facei);
                    //!!!
                    BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                    Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                    Uchemical = mixture.Hc();
                    BoundaryHE[facei] = Uthermodynamical + Uchemical;
                    //>>>
                    BoundaryT[facei] = mixture.THE(BoundaryHE[facei], BoundaryP[facei], BoundaryT[facei]);
                    Foam::Pout << "(new boundary THE U<0)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, rho = " << BoundaryRho[facei] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
                }
                // Заполнение остальных параметров
                BoundaryPsi[facei] = mixture.psi(BoundaryP[facei], BoundaryT[facei]);
                BoundaryRho[facei] = mixture.rho(BoundaryP[facei], BoundaryT[facei]);
                pmu[facei] = mixture.mu(BoundaryP[facei], BoundaryT[facei]);
                palpha[facei] = mixture.alphah(BoundaryP[facei], BoundaryT[facei]);
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::calculategetT()
{
    const Foam::fvMesh & mesh = this->p_.mesh();
    Foam::scalarField & hCells = this->he();
    const Foam::scalarField & pCells = this->p_;
    Foam::scalarField & TCells = this->T_.primitiveFieldRef();
    Foam::scalarField & psiCells = this->psi_.primitiveFieldRef();
    Foam::scalarField & rhoCells = this->rho_.primitiveFieldRef();
    Foam::scalarField & muCells = this->mu_.primitiveFieldRef();
    Foam::scalarField & alphaCells = this->alpha_.primitiveFieldRef();
    forAll(TCells, celli)
    {
        typename MixtureType::thermoType & mixture = const_cast<typename MixtureType::thermoType &>(this->cellMixture(celli));
        if(hCells[celli] - mixture.Hc() > 0.0) // Uthermodynamical
        {
            TCells[celli] = mixture.GetT(hCells[celli] - mixture.Hc(), pCells[celli]);
            if(TCells[celli] < Tmin)
            {
                const Foam::point & Position = mesh.C()[celli];
                Foam::Pout << "(old getT)T < " << Tmin << "K!, p = " << pCells[celli]/101325.0 << "atm, T = " << TCells[celli] << "K, rho = " << rhoCells[celli] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
                //<<<
                this->p_[celli] = this->Paverage[celli];
                if(this->rho_[celli] < 0.0) this->rho_[celli] = this->Rhoaverage[celli];
//if(this->rho_[celli] > RhoLimit) this->rho_[celli] = RhoLimit;
                this->T_[celli] = this->p_[celli] / this->rho_[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
                Foam::scalar Uthermodynamical = this->p_[celli] / this->rho_[celli] * mixture.Mixture->GetA(mixture.W(), this->T_[celli]);
                Foam::scalar Uchemical = hCells[celli] - Uthermodynamical;
                // Получение нового химического элемента
                std::tie(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]), Uchemical);
                mixture = this->cellMixture(celli);
                //!!!
                this->T_[celli] = this->p_[celli] / this->rho_[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
                Uthermodynamical = this->p_[celli] / this->rho_[celli] * mixture.Mixture->GetA(mixture.W(), this->T_[celli]);
                Uchemical = mixture.Hc();
                hCells[celli] = Uthermodynamical + Uchemical;
                //>>>
                TCells[celli] = mixture.GetT(hCells[celli] - mixture.Hc(), pCells[celli]);
                Foam::Pout << "(new getT)T < " << Tmin << "K!, p = " << pCells[celli]/101325.0 << "atm, T = " << TCells[celli] << "K, rho = " << rhoCells[celli] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
            }
        }else
        {
            const Foam::point & Position = mesh.C()[celli];
            Foam::Pout << "(old getT U<0)T < " << Tmin << "K!, p = " << pCells[celli]/101325.0 << "atm, T = " << TCells[celli] << "K, rho = " << rhoCells[celli] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
            //<<<
            this->p_[celli] = this->Paverage[celli];
            if(this->rho_[celli] < 0.0) this->rho_[celli] = this->Rhoaverage[celli];
//if(this->rho_[celli] > RhoLimit) this->rho_[celli] = RhoLimit;
            this->T_[celli] = this->p_[celli] / this->rho_[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
            Foam::scalar Uthermodynamical = this->p_[celli] / this->rho_[celli] * mixture.Mixture->GetA(mixture.W(), this->T_[celli]);
            Foam::scalar Uchemical = hCells[celli] - Uthermodynamical;
            // Получение нового химического элемента
            std::tie(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]), Uchemical);
            mixture = this->cellMixture(celli);
            //!!!
            this->T_[celli] = this->p_[celli] / this->rho_[celli] * mixture.W() / Foam::constant::thermodynamic::RR;
            Uthermodynamical = this->p_[celli] / this->rho_[celli] * mixture.Mixture->GetA(mixture.W(), this->T_[celli]);
            Uchemical = mixture.Hc();
            hCells[celli] = Uthermodynamical + Uchemical;
            //>>>
            TCells[celli] = mixture.GetT(hCells[celli] - mixture.Hc(), pCells[celli]);
            Foam::Pout << "(new getT U<0)T < " << Tmin << "K!, p = " << pCells[celli]/101325.0 << "atm, T = " << TCells[celli] << "K, rho = " << rhoCells[celli] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
        }
        // Заполнение остальных параметров
        psiCells[celli] = mixture.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture.rho(pCells[celli], TCells[celli]);
        muCells[celli] = mixture.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture.alphah(pCells[celli], TCells[celli]);
    }

    // Граничные условия
    Foam::volScalarField::Boundary & pBoundaryField = this->p_.boundaryFieldRef();
    Foam::volScalarField::Boundary & TBoundaryField = this->T_.boundaryFieldRef();
    Foam::volScalarField::Boundary & psiBoundaryField = this->psi_.boundaryFieldRef();
    Foam::volScalarField::Boundary & rhoBoundaryField = this->rho_.boundaryFieldRef();
    Foam::volScalarField::Boundary & heBoundaryField = this->he().boundaryFieldRef();
    Foam::volScalarField::Boundary & muBoundaryField = this->mu_.boundaryFieldRef();
    Foam::volScalarField::Boundary & alphaBoundaryField = this->alpha_.boundaryFieldRef();
    forAll(this->T_.boundaryField(), patchi)
    {
        Foam::fvPatchScalarField & BoundaryP = pBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryT = TBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryPsi = psiBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryRho = rhoBoundaryField[patchi];
        Foam::fvPatchScalarField & BoundaryHE = heBoundaryField[patchi];
        Foam::fvPatchScalarField & pmu = muBoundaryField[patchi];
        Foam::fvPatchScalarField & palpha = alphaBoundaryField[patchi];

        /*if (BoundaryT.fixesValue())
        {
            forAll(BoundaryT, facei)
            {
                const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
                BoundaryHE[facei] = mixture.HE(BoundaryP[facei], BoundaryT[facei]);
                ppsi[facei] = mixture.psi(BoundaryP[facei], BoundaryT[facei]);
                prho[facei] = mixture.rho(BoundaryP[facei], BoundaryT[facei]);
                pmu[facei] = mixture.mu(BoundaryP[facei], BoundaryT[facei]);
                palpha[facei] = mixture.alphah(BoundaryP[facei], BoundaryT[facei]);
            }
        }
        else*/
        {
            forAll(BoundaryT, facei)
            {
                typename MixtureType::thermoType & mixture = const_cast<typename MixtureType::thermoType &>(this->patchFaceMixture(patchi, facei));
                if(BoundaryHE[facei] - mixture.Hc() > 0.0) // Uthermodynamical
                {
                    BoundaryT[facei] = mixture.GetT(BoundaryHE[facei] - mixture.Hc(), BoundaryP[facei]);
                    if(BoundaryT[facei] < Tmin)
                    {
                        const Foam::point & Position = mesh.boundary()[patchi].Cf()[facei];
                        Foam::Pout << "(old boundary getT)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, mu = " << mixture.W() << ", Position = " << Position << endl;
                        //<<<
                        BoundaryP[facei] = this->Paverage.boundaryField()[patchi][facei];
                        if(BoundaryRho[facei] < 0.0) BoundaryRho[facei] = this->Rhoaverage.boundaryField()[patchi][facei];
//if(BoundaryRho[facei] > RhoLimit) BoundaryRho[facei] = RhoLimit;
                        BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                        Foam::scalar Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                        Foam::scalar Uchemical = BoundaryHE[facei] - Uthermodynamical;
                        // Получение нового химического элемента
                        std::tie(MixtureType::Induction_.boundaryFieldRef()[patchi][facei], MixtureType::MolWeight_.boundaryFieldRef()[patchi][facei]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_.boundaryField()[patchi][facei], MixtureType::MolWeight_.boundaryField()[patchi][facei]), Uchemical);
                        mixture = this->patchFaceMixture(patchi, facei);
                        //!!!
                        BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                        Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                        Uchemical = mixture.Hc();
                        BoundaryHE[facei] = Uthermodynamical + Uchemical;
                        //>>>
                        BoundaryT[facei] = mixture.GetT(Uthermodynamical, BoundaryP[facei]);
                        Foam::Pout << "(new boundary getT)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, mu = " << mixture.W() << ", Position = " << Position << endl;
                    }
                }else
                {
                    const Foam::point & Position = mesh.boundary()[patchi].Cf()[facei];
                    Foam::Pout << "(old boundary getT U<0)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, rho = " << BoundaryRho[facei] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
                    //<<<
                    BoundaryP[facei] = this->Paverage.boundaryField()[patchi][facei];
                    if(BoundaryRho[facei] < 0.0) BoundaryRho[facei] = this->Rhoaverage.boundaryField()[patchi][facei];
//if(BoundaryRho[facei] > RhoLimit) BoundaryRho[facei] = RhoLimit;
                    BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                    Foam::scalar Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                    Foam::scalar Uchemical = BoundaryHE[facei] - Uthermodynamical;
                    // Получение нового химического элемента
                    std::tie(MixtureType::Induction_.boundaryFieldRef()[patchi][facei], MixtureType::MolWeight_.boundaryFieldRef()[patchi][facei]) = GetChemistryFromHchemical(std::make_tuple(MixtureType::Induction_.boundaryField()[patchi][facei], MixtureType::MolWeight_.boundaryField()[patchi][facei]), Uchemical);
                    mixture = this->patchFaceMixture(patchi, facei);
                    //!!!
                    BoundaryT[facei] = BoundaryP[facei] / BoundaryRho[facei] * mixture.W() / Foam::constant::thermodynamic::RR;
                    Uthermodynamical = BoundaryP[facei] / BoundaryRho[facei] * mixture.Mixture->GetA(mixture.W(), BoundaryT[facei]);
                    Uchemical = mixture.Hc();
                    BoundaryHE[facei] = Uthermodynamical + Uchemical;
                    //>>>
                    BoundaryT[facei] = mixture.GetT(Uthermodynamical, BoundaryP[facei]);
                    Foam::Pout << "(new boundary getT U<0)T < " << Tmin << "K!, p = " << BoundaryP[facei]/101325.0 << "atm, T = " << BoundaryT[facei] << "K, rho = " << BoundaryRho[facei] << ", mu = " << mixture.W() << ", Position = " << Position << endl;
                }
                // Заполнение остальных параметров
                BoundaryPsi[facei] = mixture.psi(BoundaryP[facei], BoundaryT[facei]);
                BoundaryRho[facei] = mixture.rho(BoundaryP[facei], BoundaryT[facei]);
                pmu[facei] = mixture.mu(BoundaryP[facei], BoundaryT[facei]);
                palpha[facei] = mixture.alphah(BoundaryP[facei], BoundaryT[facei]);
            }
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRDEThermo<BasicPsiThermo, MixtureType>::heRDEThermo(const fvMesh& mesh, const word& phaseName):
heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    switch(UpLayer)
    {
        case EUpLayer::THE : calculateTHE();break;
        case EUpLayer::getT : calculategetT();break;
        default : return;
    };
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRDEThermo<BasicPsiThermo, MixtureType>::~heRDEThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }
    switch(UpLayer)
    {
        case EUpLayer::THE : calculateTHE();break;
        case EUpLayer::getT : calculategetT();break;
        default : return;
    };
    if (debug)
    {
        Info << "    Finished" << endl;
    }
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectRho()
{
    const Foam::scalarField & pCells = this->p_;
    const Foam::scalarField & TCells = this->T_;
    Foam::scalarField & hCells = this->he().primitiveFieldRef();
    Foam::scalarField & rhoCells = this->rho_.primitiveFieldRef();
    forAll(rhoCells, celli)
    {
        const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
        rhoCells[celli] = mixture.rho(pCells[celli], TCells[celli]);
        hCells[celli] = mixture.HE(pCells[celli], TCells[celli]);
    }
    return this->rho_;
}

template<class BasicPsiThermo, class MixtureType>
const Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Rho() const
{
    return this->rho_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::MolWeight()
{
    return MixtureType::MolWeight_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::FMolWeight()
{
    return MixtureType::FMolWeight_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Induction()
{
    return MixtureType::Induction_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::FInduction()
{
    return MixtureType::FInduction_;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::chemistry(const Foam::scalar CameraLength)
{
    const Foam::scalarField & TCells = this->T_.primitiveField();
    const Foam::scalarField & rhoCells = this->rho_.primitiveField();
    const Foam::fvMesh & mesh = this->T_.mesh();
    forAll(TCells, celli)
    {
        const Foam::point & Position = mesh.C()[celli];
        if(Position.z() < CameraLength) MixtureType::Chemistry(celli, rhoCells[celli], TCells[celli]);
        else MixtureType::Inert(celli);
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectChemistry()
{
    Foam::volScalarField & Field = MolWeight();
    forAll(Field, celli)
    {
        Foam::ChemistryElement Element = std::make_tuple(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]);
        MixtureType::CorrectChemistry(Element);
        std::tie(MixtureType::Induction_[celli], MixtureType::MolWeight_[celli]) = Element;
    }
    forAll(MixtureType::Induction_.boundaryFieldRef(), patchi) forAll (MixtureType::Induction_.boundaryFieldRef()[patchi], facei)
    {
        Foam::ChemistryElement Element = std::make_tuple(MixtureType::Induction_.boundaryField()[patchi][facei], MixtureType::MolWeight_.boundaryField()[patchi][facei]);
        MixtureType::CorrectChemistry(Element);
        std::tie(MixtureType::Induction_.boundaryFieldRef()[patchi][facei], MixtureType::MolWeight_.boundaryFieldRef()[patchi][facei]) = Element;
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectChemistry(ChemistryElement & Element)
{
    MixtureType::CorrectChemistry(Element);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Dcj()
{
    return MixtureType::CJ_D;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectErrors(bool UseChemistry, Foam::volScalarField & rho)
{
    const Foam::fvMesh & mesh = rho.mesh();
    const Foam::volVectorField & Centres = mesh.C();
    const Foam::scalarField & hCells = this->he();
    Foam::scalarField & rhoCells = this->rho_.primitiveFieldRef();
    Foam::scalarField & TCells = this->T_.primitiveFieldRef();
    Foam::volScalarField & MolWeight = this->MolWeight();
    Foam::volScalarField & Induction = this->Induction();
    Foam::scalarField & pCells = this->p_.primitiveFieldRef();
    Foam::scalarField & psiCells = this->psi_.primitiveFieldRef();
    forAll(TCells, celli)
    {
        if(TCells[celli] < Tmin || rho[celli] < 0.0)
        {
            const Foam::point & Centre = Centres[celli];
            Foam::Pout << "T < " << Tmin << "K!(celli = " << celli << ")Pold = " << pCells[celli]/101325.0 << "atm, Told = " << TCells[celli] << ", MolWeightold = " << MolWeight[celli] << ", point = " << Centre;
            //-----------------------------------------------------------
            pCells[celli] = this->Paverage[celli];
            if(rho[celli] < 0.0)
            {
                rho[celli] = this->Rhoaverage[celli];
                Foam::Pout << "(ERROR! RHO < 0) ";
                if(rho[celli] < 0.0)
                {
                    Foam::Pout << "ERROR! RHO < 0 AGAIN!" << Foam::endl;
                    rho[celli] = rho[celli+1];
                }
            }
            // установим молярную массу и температуру
            const typename MixtureType::thermoType & mixture = MixtureType::cellT(TCells[celli], celli, hCells[celli], rho[celli], pCells[celli]);
            //-------------------------------------------------------------
            psiCells[celli] = mixture.psi(pCells[celli], TCells[celli]);
            rhoCells[celli] = mixture.rho(pCells[celli], TCells[celli]);
            Foam::Pout << ", P = " << pCells[celli]/101325.0 << "atm, rho = " << rho[celli] << ", Induction = " << Induction[celli] << ", MolWeight = " << MolWeight[celli] << ", T = " << TCells[celli] << Foam::endl;
        }
    }

    // Граничные условия
    Foam::volScalarField::Boundary & pBf = this->p_.boundaryFieldRef();
    Foam::volScalarField::Boundary & TBf = this->T_.boundaryFieldRef();
    Foam::volScalarField::Boundary & MolWeightBf = MixtureType::MolWeight_.boundaryFieldRef();
    Foam::volScalarField::Boundary & InductionBf = MixtureType::Induction_.boundaryFieldRef();
    Foam::volScalarField::Boundary & psiBf = this->psi_.boundaryFieldRef();
    Foam::volScalarField::Boundary & rhoBf = this->rho_.boundaryFieldRef();
    Foam::volScalarField::Boundary & heBf = this->he().boundaryFieldRef();
    Foam::volScalarField::Boundary & muBf = this->mu_.boundaryFieldRef();
    Foam::volScalarField::Boundary & alphaBf = this->alpha_.boundaryFieldRef();
    forAll(this->T_.boundaryField(), patchi)
    {
        Foam::fvPatchScalarField & pp = pBf[patchi];
        Foam::fvPatchScalarField & pT = TBf[patchi];
        Foam::fvPatchScalarField & pMolWeight = MolWeightBf[patchi];
        Foam::fvPatchScalarField & pInduction = InductionBf[patchi];
        Foam::fvPatchScalarField & ppsi = psiBf[patchi];
        Foam::fvPatchScalarField & prho = rhoBf[patchi];
        Foam::fvPatchScalarField & phe = heBf[patchi];
        Foam::fvPatchScalarField & pmu = muBf[patchi];
        Foam::fvPatchScalarField & palpha = alphaBf[patchi];

        /*if (pT.fixesValue())
        {
            forAll(pT, facei) if(pT[facei] < Tmin || prho[facei] < 0.0)
            {
                Foam::Pout << "T < " << Tmin << "K!; BOUNDS: patchi = " << patchi << ", facei = " << facei << " Pold = " << pp[facei]/101325.0 << ", Told = " << pT[facei] << ", MolWeightold = " << pMolWeight[facei];
                //-------------------------------------------------------------
                pp[facei] = Pav.ref().boundaryFieldRef()[patchi][facei];
                // установим молярную массу и температуру
                const typename MixtureType::thermoType & mixture = MixtureType::patchFaceT(pT[facei], patchi, facei, phe[facei], prho[facei], pp[facei]);
                //-------------------------------------------------------------
                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                prho[facei] = mixture.rho(pp[facei], pT[facei]);
                Foam::Pout << ", P = " << pp[facei]/101325.0 << "atm, rho = " << prho[facei] << ", Induction = " << pInduction[facei] << ", MolWeight = " << pMolWeight[facei] << ", T = " << pT[facei] << Foam::endl;
            }
        }
        else*/
        {
            forAll(pT, facei) if(pT[facei] < Tmin || prho[facei] < 0.0)
            {
                const Foam::point & Position = mesh.boundary()[patchi].Cf()[facei];
                Foam::Pout << "T < " << Tmin << "K!; BOUNDS: patchi = " << patchi << ", facei = " << facei << " Pold = " << pp[facei]/101325.0 << ", Told = " << pT[facei] << ", MolWeightold = " << pMolWeight[facei] << ", Position = " << Position;
                //-------------------------------------------------------------
                pp[facei] = this->Paverage.boundaryFieldRef()[patchi][facei];
                // установим молярную массу и температуру
                const typename MixtureType::thermoType & mixture = MixtureType::patchFaceT(pT[facei], patchi, facei, phe[facei], prho[facei], pp[facei]);
                //-------------------------------------------------------------
                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                prho[facei] = mixture.rho(pp[facei], pT[facei]);
                Foam::Pout << ", P = " << pp[facei]/101325.0 << "atm, rho = " << prho[facei] << ", Induction = " << pInduction[facei] << ", MolWeight = " << pMolWeight[facei] << ", T = " << pT[facei] << Foam::endl;
            }
        }
    }
    //-------------------------------------------------------------
    forAll(TCells, celli) if(TCells[celli] < Tmin || pCells[celli] < 0.0 || rho[celli] < 0.0)
    {
        Foam::scalar T = TCells[celli];
        Foam::scalar P = pCells[celli];
        Foam::scalar rhorho = rho[celli];
        Foam::Pout << "ERROR! AGAIN! celli = " << celli << " T = " << T << " , P = " << P/101325.0 << " , rho = " << rhorho << Foam::endl;
    }
    // Проверка на темпетатуру (T < Tmax)
    forAll(TCells, celli)
    {
        const Foam::scalar H = this->Hthermodynamical(celli);
        Foam::scalar Tmax = this->GetTmax(MolWeight[celli], H);
        if(TCells[celli] > Tmax)
        {

            TCells[celli] = Tmax;
            const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
            Foam::scalar P = rhoCells[celli]*TCells[celli]*Foam::constant::thermodynamic::RR/MolWeight[celli];
            psiCells[celli] = mixture.psi(P, TCells[celli]);
        }
    }
}

template<class BasicPsiThermo, class MixtureType>
Foam::ChemistryElement Foam::heRDEThermo<BasicPsiThermo, MixtureType>::GetChemistry(Foam::label celli)
{
    Foam::scalar MolWeight = this->MolWeight().primitiveField()[celli];
    Foam::scalar Induction = this->Induction().primitiveField()[celli];
    return std::make_tuple(Induction, MolWeight);
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::SetChemistry(Foam::label celli, const ChemistryElement & Element)
{
    std::tie(this->Induction()[celli], this->MolWeight()[celli]) = Element;
}

template<class BasicPsiThermo, class MixtureType>
Foam::ChemistryElement Foam::heRDEThermo<BasicPsiThermo, MixtureType>::GetChemistryFromHchemical(const Foam::ChemistryElement & old, Foam::scalar Hchemical)
{
    Foam::ChemistryElement Element = MixtureType::GetFromHchemical(old, Hchemical);
    MixtureType::CorrectChemistry(Element);
    return Element;
}

template<class BasicPsiThermo, class MixtureType>
Foam::ChemistryElement Foam::heRDEThermo<BasicPsiThermo, MixtureType>::GetChemistryFromHpotential(Foam::scalar rho, Foam::scalar P, Foam::scalar gamma, Foam::scalar Hpotential, const Foam::ChemistryElement & old)
{
    Foam::scalar Induction = std::get<0>(old);
    Foam::scalar MolWeight = std::get<1>(old);
    const scalar Hchemical = Hpotential - P/rho/(gamma - 1.0);
    Foam::ChemistryElement Element = MixtureType::GetFromHchemical(old, Hchemical);
    return Element;
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Hthermodynamical(Foam::scalar P, Foam::scalar rho, const Foam::ChemistryElement & Element)
{
    Foam::scalar Induction = std::get<0>(Element);
    Foam::scalar MolWeight = std::get<1>(Element);
    Foam::scalar T = P / rho * MolWeight / Foam::constant::thermodynamic::RR;
    Foam::scalar Hthermodinamical = P/rho*MixtureType::GetA(MolWeight, T);
    return Hthermodinamical;
    //const typename MixtureType::thermoType & mixture = MixtureType::Mixture(Element);
    //return mixture.Hs(P,T);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Hchemical(const Foam::ChemistryElement & Element)
{
    const typename MixtureType::thermoType & mixture = MixtureType::Mixture(Element);
    return mixture.Hc();
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Hpotential(Foam::scalar P, Foam::scalar rho, const Foam::ChemistryElement & Element)
{
    Foam::scalar Induction = std::get<0>(Element);
    Foam::scalar MolWeight = std::get<1>(Element);
    Foam::scalar T = P / rho * MolWeight / Foam::constant::thermodynamic::RR;
    const typename MixtureType::thermoType & mixture = MixtureType::Mixture(Element);
    Foam::scalar Hchemical = mixture.Hc();
    Foam::scalar Hthermodinamical = P/rho*MixtureType::GetA(MolWeight, T);
    return Hchemical + Hthermodinamical;
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Hthermodynamical(Foam::label celli)
{
    const Foam::scalarField & TCells = this->T_.primitiveField();
    const Foam::scalarField & pCells = this->p_.primitiveField();
    const Foam::scalarField & hCells = this->he();
    const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
    return hCells[celli] - mixture.Hc();
    //return mixture.Hs(pCells[celli], TCells[celli]);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Hchemical(Foam::label celli)
{
    const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
    return mixture.Hc();
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::MolWeightStar()
{
    return MixtureType::MolWeight0;
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::patchInletButtMolWeight(const label patchi, const label facei)
{
    const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
    return MixtureType::InletButtMolWeight(patchi, facei, mixture);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::patchInletSideWallMolWeight(const label patchi, const label facei)
{
    const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
    return MixtureType::InletSideWallMolWeight(patchi, facei, mixture);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::patchOutletMolWeight(const label patchi, const label facei)
{
    const Foam::volScalarField::Boundary & BoundaryMolWeight = MixtureType::MolWeight_.boundaryField();
    return BoundaryMolWeight[patchi][facei];
}

template<class BasicPsiThermo, class MixtureType>
const Foam::word Foam::heRDEThermo<BasicPsiThermo, MixtureType>::InfoElement(const label celli)
{
    return MixtureType::InfoElement(celli);
}

template<class BasicPsiThermo, class MixtureType>
const Foam::word Foam::heRDEThermo<BasicPsiThermo, MixtureType>::InfoElement(const label patchi, const label facei)
{
    return MixtureType::InfoElement(patchi, facei);
}

#endif
