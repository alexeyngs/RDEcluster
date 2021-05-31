
#ifndef MolWeightInletButtBoundary_H
#define MolWeightInletButtBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
    
class MolWeightInletButtBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("MolWeightInletButtBoundary");
    // Constructors
    //- Construct from patch and internal field
    MolWeightInletButtBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    MolWeightInletButtBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    MolWeightInletButtBoundary(const MolWeightInletButtBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    MolWeightInletButtBoundary(const MolWeightInletButtBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new MolWeightInletButtBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    MolWeightInletButtBoundary(const MolWeightInletButtBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new MolWeightInletButtBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~MolWeightInletButtBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
