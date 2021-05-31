#include "TInletButtBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "thermodynamicConstants.H"
#include "rhoReactionThermo.H"
#include "../solver/thermo/rhoRDEThermo.H"
namespace Foam
{
// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
//{{{ begin localCode
//}}} end localCode
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void TInletButtBoundary_1_0(bool load)
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
    TInletButtBoundary
);
const char * const TInletButtBoundary::SHA1sum = "66d4e0b38a18f2b7674332f6d021721a92f0e80e";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
TInletButtBoundary::TInletButtBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF):fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct TInletButtBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField\n";
    }
}

TInletButtBoundary::TInletButtBoundary(const TInletButtBoundary& ptf,const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const fvPatchFieldMapper& mapper):fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct TInletButtBoundary from patch/DimensionedField/mapper" << endl;
    }
}

TInletButtBoundary::TInletButtBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const dictionary& dict):fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct TInletButtBoundary from patch/dictionary" << endl;
    }
}

TInletButtBoundary::TInletButtBoundary(const TInletButtBoundary& ptf):fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct TInletButtBoundary as copy" << endl;;
    }
}

TInletButtBoundary::TInletButtBoundary(const TInletButtBoundary& ptf, const DimensionedField<scalar, volMesh>& iF):fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct TInletButtBoundary as copy/DimensionedField\n";
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
TInletButtBoundary::~TInletButtBoundary()
{
    if (false)
    {
        Info << "destroy TInletButtBoundary" << endl;
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void TInletButtBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info << "updateCoeffs TInletButtBoundary" << endl;
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
        const Foam::scalar MolWeight = Thermo.patchInletButtMolWeight(ThisPatch, facei);
        Foam::scalar MolWeightStar = Thermo.MolWeightStar();
        result[facei] = ::TinButt(scalarParameters, time, MolWeightStar, MolWeight, p.ref()[facei], T.ref()[facei], rho.ref()[facei], U.ref()[facei], c.ref()[facei], gamma.ref()[facei]);
    }
//this->operator==(300.0);
// end code
    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}
} // End namespace Foam


/*
    Пример граничной ячейки:
    const label patchi = this->patch().template index();
    const scalarField TempIntF(Temp.boundaryField()[patchi].patchInternalField());
    const volScalarField Temp = db().lookupObject<volScalarField>(TempName_);
    или так:
    const label patchi = this->patch().template index();
    const volScalarField Temp = this->db().template lookupObject<volScalarField>(TempName_);
    Пример на границе:
    const scalarField & rho = this->patch().lookupPatchField<volScalarField, scalar>("rho");
*/