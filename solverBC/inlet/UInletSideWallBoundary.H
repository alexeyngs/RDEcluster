#ifndef UInletSideWallBoundary_H
#define UInletSideWallBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class UInletSideWallBoundary : public fixedValueFvPatchField<vector>
{
public:
    //- Information about the SHA1 of the code itself
    static const char * const SHA1sum;

    //- Runtime type information
    TypeName("UInletSideWallBoundary");

    // Constructors
    //- Construct from patch and internal field
    UInletSideWallBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&);

    //- Construct from patch, internal field and dictionary
    UInletSideWallBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&,const dictionary&);

    //- Construct by mapping a copy onto a new patch
    UInletSideWallBoundary(const UInletSideWallBoundary&,const fvPatch&,const DimensionedField<vector, volMesh>&,const fvPatchFieldMapper&);

    //- Construct as copy
    UInletSideWallBoundary(const UInletSideWallBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField > clone() const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletSideWallBoundary(*this)
        );
    }

    //- Construct as copy setting internal field reference
    UInletSideWallBoundary(const UInletSideWallBoundary&,const DimensionedField<vector, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletSideWallBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~UInletSideWallBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif