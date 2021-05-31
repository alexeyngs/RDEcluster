
#ifndef RDEspecie_H
#define RDEspecie_H
#include "word.H"
#include "scalar.H"
#include "dictionary.H"
#include "thermodynamicConstants.H"

namespace Foam
{
// Forward declaration of friend functions and operators
class RDEspecie;
inline RDEspecie operator + (const RDEspecie&, const RDEspecie&);
inline RDEspecie operator * (const scalar, const RDEspecie&);
inline RDEspecie operator == (const RDEspecie&, const RDEspecie&);
Ostream & operator << (Ostream&, const RDEspecie&);

class RDEspecie
{
public:

        //- Name of specie
        word name_;

        //- Number of moles of this component in the mixture
        scalar moles_;

        // - Количество атомов в смеси
        scalar atoms_;

        //- Количество атомов, на которых будет разлагаться смесь
        scalar decomposition_;

        //- Количество атомов кислорода, которое необходимо для полного сжигания компонента
        scalar needO_;

        //- Молекулярная масса [kg/kmol]
        mutable scalar MolWeight_;

        //- Энергия диссоциации компонента смеси        
        mutable scalar Edissociation_;

        //- Энергия активации [Дж/моль]
        mutable scalar Einduction_;

        //- Константа активации в моль*с/м3
        mutable scalar Kinduction_;

        //- Показатель степени в концентрации
        mutable scalar InductionPower_;

        //- Параметр тета - характерная температура возбуждения колебательных степеней свободы
        mutable scalar theta_;
public:
    //- Runtime type information
    ClassName("RDEspecie");
    // Constructors
    //- Construct from components without name
    inline RDEspecie(const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight, const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction, const scalar InductionPower, const scalar theta);

    //- Construct from components with name
    inline RDEspecie(const word & name, const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight, const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction, const scalar InductionPower, const scalar theta);

    //- Construct as named copy
    inline RDEspecie(const word & name, const RDEspecie&);

    //- Construct from dictionary
    RDEspecie(const dictionary & dict);
// Member Functions
    //- Name
    inline const word & name() const;

    //- Молекулярная масса [kg/kmol]
    inline scalar W() const;

    //- Количество атомов кислорода, которое необходимо для полного сжигания компонента
    inline scalar needO() const;

    //- Количество молей
    inline scalar moles() const;

    //- Количество атомов в смеси
    inline scalar atoms() const;

    //- Количество атомов, на которых будет разлагаться смесь
    inline scalar decomposition() const;

    //- Энергия диссоциации компонента смеси        
    inline scalar Edissociation() const;

    //- Энергия активации [Дж/моль]
    inline scalar Einduction() const;

    //- Константа активации в моль*с/м3
    inline scalar Kinduction() const;

    //- Показатель степени в концентрации
    inline scalar InductionPower() const;

    //- температура возбуждения колебательных степеней свободы
    inline scalar theta() const;

    //- No of moles of this RDEspecies in mixture
    inline scalar Y() const;

    //- Gas constant [J/(kg K)]
    inline scalar R() const;

//- Write to Ostream
    void write(Ostream& os) const;
    // Member operators
    inline void operator = (const RDEspecie &);
    inline void operator += (const RDEspecie &);
    inline void operator *= (const scalar);

    // Friend operators
    inline friend RDEspecie operator + (const RDEspecie &, const RDEspecie &);
    inline friend RDEspecie operator * (const scalar, const RDEspecie &);
    inline friend RDEspecie operator == (const RDEspecie &, const RDEspecie &);
    // Ostream Operator
    friend Ostream & operator << (Ostream &, const RDEspecie &);
};
} // End namespace Foam

// Конструкторы
inline Foam::RDEspecie::RDEspecie(const Foam::word & name, const Foam::scalar moles, const Foam::scalar atoms, const Foam::scalar decomposition, const Foam::scalar molWeight,
    const Foam::scalar needO, const Foam::scalar Edissociation, const Foam::scalar Einduction, const Foam::scalar Kinduction, const scalar InductionPower, const scalar theta):
    name_(name), moles_(moles), atoms_(atoms), decomposition_(decomposition), MolWeight_(molWeight), needO_(needO),
    Edissociation_(Edissociation), Einduction_(Einduction), Kinduction_(Kinduction), InductionPower_(InductionPower), theta_(theta)
{}

inline Foam::RDEspecie::RDEspecie(const Foam::scalar moles, const Foam::scalar atoms, const Foam::scalar decomposition, const Foam::scalar molWeight,
    const Foam::scalar needO, const Foam::scalar Edissociation, const Foam::scalar Einduction, const Foam::scalar Kinduction, const scalar InductionPower, const scalar theta):
    moles_(moles), atoms_(atoms), decomposition_(decomposition), MolWeight_(molWeight), needO_(needO),
    Edissociation_(Edissociation), Einduction_(Einduction), Kinduction_(Kinduction), InductionPower_(InductionPower), theta_(theta)
{}

inline Foam::RDEspecie::RDEspecie(const Foam::word & name, const Foam::RDEspecie & st):
    name_(name), moles_(st.moles_), MolWeight_(st.MolWeight_), needO_(st.needO_), atoms_(st.atoms_),
    Edissociation_(st.Edissociation_), Einduction_(st.Einduction_), Kinduction_(st.Kinduction_), InductionPower_(st.InductionPower_), theta_(st.theta_),
    decomposition_(st.decomposition_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word & Foam::RDEspecie::name() const
{
    return name_;
}

inline Foam::scalar Foam::RDEspecie::W() const
{
    return MolWeight_;
}

inline Foam::scalar Foam::RDEspecie::moles() const
{
    return moles_;
}

inline Foam::scalar Foam::RDEspecie::atoms() const
{
    return atoms_;
}

inline Foam::scalar Foam::RDEspecie::decomposition() const
{
    return decomposition_;
}

inline Foam::scalar Foam::RDEspecie::Edissociation() const
{
    return Edissociation_;
}

inline Foam::scalar Foam::RDEspecie::Einduction() const
{
    return Einduction_;
}

inline Foam::scalar Foam::RDEspecie::Kinduction() const
{
    return Kinduction_;
}

inline Foam::scalar Foam::RDEspecie::InductionPower() const
{
    return InductionPower_;
}

inline Foam::scalar Foam::RDEspecie::theta() const
{
    return theta_;
}

inline Foam::scalar Foam::RDEspecie::needO() const
{
    return needO_;
}

inline Foam::scalar Foam::RDEspecie::Y() const
{
    return 1.0;
}


inline Foam::scalar Foam::RDEspecie::R() const
{
    return Foam::constant::thermodynamic::RR/MolWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::RDEspecie::operator = (const Foam::RDEspecie & st)
{
    //name_ = st.name_;
    moles_ = st.moles_;
    atoms_ = st.atoms_;
    decomposition_ = st.decomposition_;
    MolWeight_ = st.MolWeight_;
    needO_ = st.needO_;
    Edissociation_ = st.Edissociation_;
    Einduction_ = st.Einduction_;
    Kinduction_ = st.Kinduction_;
    theta_ = st.theta_;
}


inline void Foam::RDEspecie::operator += (const Foam::RDEspecie & st)
{
    const scalar summoles = moles_ + st.moles_;
    // st.moles_ - сумма старых молярных масс
    // summoles - сумма новых молярных масс
    if (mag(summoles) > SMALL)
    {
        atoms_ = (moles_*atoms_ + st.moles_*st.atoms_) / summoles;
        decomposition_ = (moles_*decomposition_ + st.moles_*st.decomposition_) / summoles;
        MolWeight_ = (moles_*MolWeight_ + st.moles_*st.MolWeight_) / summoles;
        needO_ = (moles_*needO_ + st.moles_*st.needO_) / summoles;
        Edissociation_ = (moles_*Edissociation_ + st.moles_*st.Edissociation_) / summoles;
        Einduction_ = (moles_*Einduction_ + st.moles_*st.Einduction_) / summoles;
        Kinduction_ = (moles_*Kinduction_ + st.moles_*st.Kinduction_) / summoles;
        theta_ = (moles_*theta_ + st.moles_*st.theta_) / summoles;
    }
    moles_ = summoles;
}

inline void Foam::RDEspecie::operator *= (const Foam::scalar s)
{
    moles_ *= s;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::RDEspecie operator + (const Foam::RDEspecie & st1, const Foam::RDEspecie & st2)
{
    const Foam::scalar summoles = Foam::max(st1.moles_ + st2.moles_, Foam::SMALL);
    if (Foam::mag(summoles) > Foam::SMALL)
    {
        return Foam::RDEspecie
        (
            summoles,
            (st1.moles_*st1.atoms_ + st2.moles_*st2.atoms_) / summoles,
            (st1.moles_*st1.decomposition_ + st2.moles_*st2.decomposition_) / summoles,
            (st1.moles_*st1.MolWeight_ + st2.moles_*st2.MolWeight_) / summoles,
            (st1.moles_*st1.needO_ + st2.moles_*st2.needO_) / summoles,
            (st1.moles_*st1.Edissociation_ + st2.moles_*st2.Edissociation_) / summoles,
            (st1.moles_*st1.Einduction_ + st2.moles_*st2.Einduction_) / summoles,
            (st1.moles_*st1.Kinduction_ + st2.moles_*st2.Kinduction_) / summoles,
            (st1.moles_*st1.theta_ + st2.moles_*st2.theta_) / summoles,
            0.0 // Показамтель степени в концентрации
        );
    }
    else
    {
        return st1;
    }
}

namespace Foam{
// Зачем-то сюда заходит...
inline RDEspecie operator * (const scalar s, const RDEspecie & st)
{
    return Foam::RDEspecie
    (
        st.moles_,
        st.atoms_,
        st.decomposition_,
        st.MolWeight_,
        st.needO_,
        st.Edissociation_,
        st.Einduction_,
        st.Kinduction_,
        st.InductionPower_,
        st.theta_
    );
}
}// end of namespace Foam

inline Foam::RDEspecie operator == (const Foam::RDEspecie & st1, const Foam::RDEspecie & st2)
{
    Foam::scalar diff_moles = st2.moles_ - st1.moles_;
    if(Foam::mag(diff_moles) < Foam::SMALL)
    {
        diff_moles = Foam::SMALL;
    }

    const Foam::scalar diff_atoms = st2.moles_*st2.atoms_ - st1.moles_*st1.atoms_;
    const Foam::scalar diff_decomposition = st2.moles_*st2.decomposition_ - st1.moles_*st1.decomposition_;
    const Foam::scalar diff_molWeight = st2.moles_*st2.MolWeight_ - st1.moles_*st1.MolWeight_;
    const Foam::scalar diff_NeedO = st2.moles_*st2.needO_ - st1.moles_*st1.needO_;
    const Foam::scalar diff_Edissociation = st2.moles_*st2.Edissociation_ - st1.moles_*st1.Edissociation_;
    const Foam::scalar diff_Einduction = st2.moles_*st2.Einduction_ - st1.moles_*st1.Einduction_;
    const Foam::scalar diff_Kinduction = st2.moles_*st2.Kinduction_ - st1.moles_*st1.Kinduction_;
    const Foam::scalar diff_InductionPower = st2.moles_*st2.InductionPower_ - st1.moles_*st1.InductionPower_;
    const Foam::scalar diff_theta = st2.moles_*st2.theta_ - st1.moles_*st1.theta_;


    
    
    Foam::scalar atoms = 0.0;
    Foam::scalar decomposition = 0.0;
    Foam::scalar molWeight = 0.0;
    Foam::scalar needO = 0.0;
    Foam::scalar Edissociation = 0.0;
    Foam::scalar Einduction = 0.0;
    Foam::scalar Kinduction = 0.0;
    Foam::scalar InductionPower = 0.0;
    Foam::scalar theta = 0.0;
    

    if (Foam::mag(diff_atoms) > Foam::SMALL) atoms = diff_atoms / diff_moles;
    if (Foam::mag(diff_decomposition) > Foam::SMALL) decomposition = diff_decomposition / diff_moles;
    if (Foam::mag(diff_molWeight) > Foam::SMALL) molWeight = diff_molWeight / diff_moles;
    if (Foam::mag(diff_NeedO) > Foam::SMALL) needO = diff_NeedO / diff_moles;
    if (Foam::mag(diff_Edissociation) > Foam::SMALL) Edissociation = diff_Edissociation / diff_moles;
    if (Foam::mag(diff_Einduction) > Foam::SMALL) Einduction = diff_Einduction / diff_moles;
    if (Foam::mag(diff_Kinduction) > Foam::SMALL) Kinduction = diff_Kinduction / diff_moles;
    if (Foam::mag(diff_InductionPower) > Foam::SMALL) InductionPower = diff_InductionPower / diff_moles;
    if (Foam::mag(diff_theta) > Foam::SMALL) theta = diff_theta / diff_moles;

    return Foam::RDEspecie(diff_moles, atoms, decomposition, molWeight, needO, Edissociation, Einduction, Kinduction, InductionPower, theta);
}

#endif
