
#ifndef RHORDETHERMO_H
#define RHORDETHERMO_H
#include "moleFractions.H"
#include "rhoReactionThermo.H"
#include <tuple>
namespace Foam
{
    typedef std::tuple<scalar, scalar> ChemistryElement;
class rhoRDEThermo : public rhoReactionThermo
{
public:
    //- Runtime type information
    TypeName("rhoRDEThermo");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoRDEThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );
    //- Construct from mesh and phase name
    rhoRDEThermo(const fvMesh&,const word& phaseName);
    //- Standard selection based on fvMesh
    static autoPtr<rhoRDEThermo> New(const fvMesh&,const word& phaseName=word::null);
    //- Destructor
    virtual ~rhoRDEThermo();
    //---------------------------------------------
    volScalarField Paverage;
    volScalarField Rhoaverage;
    //---------------------------------------------
    virtual volScalarField & CorrectRho() = 0;
    virtual volScalarField & MolWeight() = 0;
    virtual volScalarField & FMolWeight() = 0;
    virtual volScalarField & Induction() = 0;
    virtual volScalarField & FInduction() = 0;
    virtual void chemistry(const scalar CameraLength) = 0;
    virtual void CorrectChemistry() = 0;
    virtual void CorrectChemistry(ChemistryElement & Element) = 0;
    virtual scalar Dcj() = 0;
    virtual void CorrectErrors(bool UseChemistry, volScalarField & rho, surfaceScalarField & Rhof, surfaceScalarField & Pf) = 0;
    virtual ChemistryElement GetChemistry(label celli) = 0;
    virtual void SetChemistry(label celli, const ChemistryElement & Element) = 0;
    virtual ChemistryElement GetChemistryFromHchemical(const ChemistryElement & old, scalar Hchemical) = 0;
    virtual ChemistryElement GetChemistryFromHpotential(scalar rho, scalar P, scalar gamma, scalar Hpotential, const ChemistryElement & old) = 0;
    virtual scalar Hthermodynamical(label celli) = 0;
    virtual scalar Hchemical(label celli) = 0;
    virtual scalar Hchemical(const ChemistryElement & Element) = 0;
    virtual scalar Hthermodynamical(scalar P, scalar rho, const ChemistryElement & Element) = 0;
    virtual scalar Hpotential(scalar P, scalar rho, const ChemistryElement & Element) = 0;
    virtual scalar MolWeightStar() = 0;
    virtual scalar patchInletButtMolWeight(const label patchi, const label facei) = 0;
    virtual scalar patchInletSideWallMolWeight(const label patchi, const label facei) = 0;
    virtual scalar patchOutletMolWeight(const label patchi, const label facei) = 0;
    virtual const word InfoElement(const label celli) = 0;
    virtual const word InfoElement(const label patchi, const label facei) = 0;
};

} // End namespace Foam

#endif // RHORDETHERMO_H
