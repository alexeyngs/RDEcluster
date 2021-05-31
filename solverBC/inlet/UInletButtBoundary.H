#ifndef UInletButtBoundary_H
#define UInletButtBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class UInletButtBoundary : public fixedValueFvPatchField<vector>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("UInletButtBoundary");

    // Constructors
    //- Construct from patch and internal field
    UInletButtBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&);

    //- Construct from patch, internal field and dictionary
    UInletButtBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&,const dictionary&);

    //- Construct by mapping a copy onto a new patch
    UInletButtBoundary(const UInletButtBoundary&,const fvPatch&,const DimensionedField<vector, volMesh>&,const fvPatchFieldMapper&);

    //- Construct as copy
    UInletButtBoundary(const UInletButtBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField > clone() const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletButtBoundary(*this)
        );
    }

    //- Construct as copy setting internal field reference
    UInletButtBoundary(const UInletButtBoundary&,const DimensionedField<vector, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletButtBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~UInletButtBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif