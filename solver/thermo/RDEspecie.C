#include "RDEspecie.H"
#include "constants.H"

namespace Foam
{
    defineTypeNameAndDebug(RDEspecie, 0);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RDEspecie::RDEspecie(const dictionary & dict):
    name_(dict.dictName()),
    moles_(readScalar(dict.subDict("specie").lookup("moles"))),
    atoms_(readScalar(dict.subDict("specie").lookup("atoms"))),
    decomposition_(readScalar(dict.subDict("specie").lookup("decomposition"))),
    MolWeight_(readScalar(dict.subDict("specie").lookup("molWeight"))),
    needO_(dict.subDict("specie").lookupOrDefault("needO", 0.0)),
    Edissociation_(dict.subDict("specie").lookupOrDefault("Edissociation", 0.0)),
    Einduction_(dict.subDict("specie").lookupOrDefault("Einduction", 0.0)),
    Kinduction_(dict.subDict("specie").lookupOrDefault("Kinduction", 0.0)),
    InductionPower_(dict.subDict("specie").lookupOrDefault("InductionPower", 0.0)),
    theta_(dict.subDict("specie").lookupOrDefault("theta", 0.0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::RDEspecie::write(Ostream& os) const
{
    dictionary dict("specie");
    dict.add("moles", moles_);
    dict.add("atoms", atoms_);
    dict.add("decomposition", decomposition_);
    dict.add("molWeight", MolWeight_);
    dict.add("needO", needO_);
    dict.add("Edissociation", Edissociation_);
    dict.add("Einduction", Einduction_);
    dict.add("Kinduction", Kinduction_);
    dict.add("InductionPower", InductionPower_);
    dict.add("theta", theta_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //
Foam::Ostream & Foam::operator << (Foam::Ostream & os, const Foam::RDEspecie & st)
{
    st.write(os);
    os.check("Ostream& operator<<(Ostream& os, const RDEspecie& st)");
    return os;
}
