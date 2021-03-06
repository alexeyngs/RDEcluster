
#ifndef MolWeightInletSideWallBoundary_H
#define MolWeightInletSideWallBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
    
class MolWeightInletSideWallBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("MolWeightInletSideWallBoundary");
    // Constructors
    //- Construct from patch and internal field
    MolWeightInletSideWallBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    MolWeightInletSideWallBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    MolWeightInletSideWallBoundary(const MolWeightInletSideWallBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    MolWeightInletSideWallBoundary(const MolWeightInletSideWallBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new MolWeightInletSideWallBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    MolWeightInletSideWallBoundary(const MolWeightInletSideWallBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new MolWeightInletSideWallBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~MolWeightInletSideWallBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
