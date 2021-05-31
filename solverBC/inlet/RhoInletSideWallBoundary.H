
#ifndef RhoInletSideWallBoundary_H
#define RhoInletSideWallBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
    
class RhoInletSideWallBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("RhoInletSideWallBoundary");
    // Constructors
    //- Construct from patch and internal field
    RhoInletSideWallBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    RhoInletSideWallBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    RhoInletSideWallBoundary(const RhoInletSideWallBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    RhoInletSideWallBoundary(const RhoInletSideWallBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new RhoInletSideWallBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    RhoInletSideWallBoundary(const RhoInletSideWallBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new RhoInletSideWallBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~RhoInletSideWallBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
