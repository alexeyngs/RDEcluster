#include "RhoInletButtBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "thermodynamicConstants.H"
#include "rhoReactionThermo.H"

#include "../solver/thermo/rhoRDEThermo.H"
#include "boundary.H"
namespace Foam
{
// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
//{{{ begin localCode
//}}} end localCode
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void RhoInletButtBoundary_1_0(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    RhoInletButtBoundary
);
const char * const RhoInletButtBoundary::SHA1sum = "66d4e0b38a18f2b7674332f6d021721a92f0e80e";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
RhoInletButtBoundary::RhoInletButtBoundary(const fvPatch & p, const DimensionedField<scalar, volMesh> & iF):
fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField" << endl;
    }
}
RhoInletButtBoundary::RhoInletButtBoundary(const RhoInletButtBoundary& ptf,const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const fvPatchFieldMapper& mapper):
fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField/mapper\n";
    }
}
RhoInletButtBoundary::RhoInletButtBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const dictionary& dict):
fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info <<"construct RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/dictionary\n";
    }
}

RhoInletButtBoundary::RhoInletButtBoundary(const RhoInletButtBoundary& ptf):
fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info <<"construct RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e as copy\n";
    }
}

RhoInletButtBoundary::RhoInletButtBoundary(const RhoInletButtBoundary& ptf,const DimensionedField<scalar, volMesh>& iF):
fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info <<"construct RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e as copy/DimensionedField\n";
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
RhoInletButtBoundary::~RhoInletButtBoundary()
{
    if (false)
    {
        Info <<"destroy RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e\n";
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void RhoInletButtBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info<<"updateCoeffs RhoInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e\n";
    }
// begin code
    const scalar & time = this->db().time().value();   // Время
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    const Foam::vectorField & Cf = this->patch().Cf();
    // Получим текущий патч - границу
    const label ThisPatch = patch().index();
    scalarField & result = *this;
    // на граничной ячейки
    tmp<scalarField> Field(this->patchInternalField());
    tmp<scalarField> rho = this->db().template lookupObject<volScalarField>("rho").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> p = this->db().template lookupObject<volScalarField>("p").boundaryField()[ThisPatch].patchInternalField();
    tmp<vectorField> U = this->db().template lookupObject<volVectorField>("U").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> c = this->db().template lookupObject<volScalarField>("c").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> gamma = this->db().template lookupObject<volScalarField>("gamma").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> T = this->db().template lookupObject<volScalarField>("T").boundaryField()[ThisPatch].patchInternalField();
    //-----------------------------------------------------
    const rhoReactionThermo & thermo = this->db().lookupObject<rhoReactionThermo>("thermophysicalProperties");
    rhoRDEThermo & Thermo = const_cast<rhoRDEThermo &> (dynamic_cast<const rhoRDEThermo &> (thermo));
    forAll(Cf, facei)
    {
        Foam::scalar MolWeight = Thermo.patchInletButtMolWeight(ThisPatch, facei);
        Foam::scalar MolWeightStar = Thermo.MolWeightStar();
        result[facei] = ::RhoinButt(scalarParameters, time, MolWeightStar, MolWeight, p.ref()[facei], T.ref()[facei], rho.ref()[facei], U.ref()[facei], c.ref()[facei], gamma.ref()[facei]);
    }
//this->operator==(101325.0);
// end code
    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}
} // End namespace Foam
