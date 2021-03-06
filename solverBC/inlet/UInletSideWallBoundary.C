#include "UInletSideWallBoundary.H"
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
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void UInletSideWallBoundary_1_0(bool load)
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

makeRemovablePatchTypeField(fvPatchVectorField,UInletSideWallBoundary);
const char* const UInletSideWallBoundary::SHA1sum = "2cb7a0c4345a61b88d329a26e832e625c53cdd92";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
UInletSideWallBoundary::UInletSideWallBoundary(const fvPatch& p,const DimensionedField<vector, volMesh>& iF):fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92"
            " from patch/DimensionedField\n";
    }
}

UInletSideWallBoundary::UInletSideWallBoundary(const UInletSideWallBoundary& ptf,const fvPatch& p,const DimensionedField<vector, volMesh>& iF,const fvPatchFieldMapper& mapper):fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct from patch/DimensionedField/mapper\n";
    }
}

UInletSideWallBoundary::UInletSideWallBoundary(const fvPatch& p,const DimensionedField<vector, volMesh>& iF, const dictionary& dict):fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92 from patch/dictionary\n";
    }
}

UInletSideWallBoundary::UInletSideWallBoundary(const UInletSideWallBoundary& ptf):fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92"
            " as copy\n";
    }
}

UInletSideWallBoundary::UInletSideWallBoundary(const UInletSideWallBoundary& ptf,const DimensionedField<vector, volMesh>& iF):fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92 "
            "as copy/DimensionedField\n";
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
UInletSideWallBoundary::~UInletSideWallBoundary()
{
    if (false)
    {
        Info<<"destroy ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92\n";
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void UInletSideWallBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info<<"updateCoeffs ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92\n";
    }
//begin code
    const scalar & time = this->db().time().value();   // Время
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    const Foam::vectorField & Cf = this->patch().Cf();
    Foam::tmp<Foam::vectorField> Normals = this->patch().nf(); // Массив нормалей
     // Получим текущий патч - границу
    const label ThisPatch = patch().index();
    vectorField & result = *this;
    // на граничной ячейки    
    tmp<vectorField> Field(this->patchInternalField());
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
        Foam::scalar MolWeight = Thermo.patchInletSideWallMolWeight(ThisPatch, facei);
        Foam::scalar MolWeightStar = Thermo.MolWeightStar();
        result[facei] = ::UinSideWall(scalarParameters, time, Normals.ref()[facei], MolWeightStar, MolWeight, p.ref()[facei], T.ref()[facei], rho.ref()[facei], U.ref()[facei], c.ref()[facei], gamma.ref()[facei]);
    }
//this->operator==(vector(0,0,0));
// end code
    this->fixedValueFvPatchField<vector>::updateCoeffs();
}

} // End namespace Foam
