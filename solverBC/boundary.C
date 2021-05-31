#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>
#include "fixedValueFvPatchFields.H"
#include "thermodynamicConstants.H"
#include "fvCFD.H"
#include "scalarIOList.H"
#include "Swap.H"
#include "../solver/interface.H"

EBreakDown TypeBreakDown = EBreakDown::None;
void ReadBreakDown(const Foam::dictionary & dict)
{
    word BreakDownString = dict.lookup("BreakDown");
    if(BreakDownString == "Accustic") TypeBreakDown = EBreakDown::Accustic;
    else if(BreakDownString == "Godunov") TypeBreakDown = EBreakDown::Godunov;
    else if(BreakDownString == "Kolgan") TypeBreakDown = EBreakDown::Kolgan;
}
//====================================================================================================================================

// Входящие параметры
Foam::scalar PinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar,  const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[5];
    const Foam::scalar & Pstar = scalarParameters[6];
    const Foam::scalar & Tstar = scalarParameters[7];
    //------------------------------------------------------------------------
    return Pcenter;
}

Foam::vector UinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::vector Normal,  const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter,  const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[5];
    const Foam::scalar & Pstar = scalarParameters[6];
    const Foam::scalar & Tstar = scalarParameters[7];
    //------------------------------------------------------------------------
    Foam::vector result;
    if(Pcenter >= Pstar)
    {
        result = Ucenter - (Ucenter & Normal) * Normal; // Тангенсальная компонента скорости
    }else
    {
        Foam::scalar F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        Foam::scalar X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        Foam::scalar U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * MolWeightStar / Foam::constant::thermodynamic::RR / Tstar;
        U *= Foam::sqrt(Pstar/Rhostar); // Выход из нормировки на *
        result = Normal*U;
    }
    return result;
}

Foam::scalar RhoinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[5];
    const Foam::scalar & Pstar = scalarParameters[6];
    const Foam::scalar & Tstar = scalarParameters[7];
    //------------------------------------------------------------------------
    Foam::scalar rho;
    Foam::scalar F, X, U;
    if(Pcenter >= Pstar)
    {
        rho = Rhocenter;
    }else
    {
        F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * MolWeightStar / Foam::constant::thermodynamic::RR / Tstar;
        rho = -F/U;
        rho *= Rhostar; // Выход из нормировки на *
    }
    return rho;
}

Foam::scalar TinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter,  const Foam::scalar Tcenter,const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[5];
    const Foam::scalar & Pstar = scalarParameters[6];
    const Foam::scalar & Tstar = scalarParameters[7];
    //------------------------------------------------------------------------
    Foam::scalar Rho = RhoinButt(scalarParameters, time, MolWeightStar, MolWeight, Pcenter, Tcenter, Rhocenter, Ucenter, c, gamma);
    Foam::scalar P = PinButt(scalarParameters, time, MolWeightStar, MolWeight, Pcenter, Tcenter, Rhocenter, Ucenter, c, gamma);
    return P * MolWeight / Rho / Foam::constant::thermodynamic::RR;
}

//------------------------------------------------------------------------
Foam::scalar PinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[8];
    const Foam::scalar & Pstar = scalarParameters[9];
    const Foam::scalar & Tstar = scalarParameters[10];
    //------------------------------------------------------------------------
    return Pcenter;
}

Foam::vector UinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::vector Normal, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[8];
    const Foam::scalar & Pstar = scalarParameters[9];
    const Foam::scalar & Tstar = scalarParameters[10];
    //------------------------------------------------------------------------
    Foam::scalar F, X, U;
    Foam::vector result;
    if(Pcenter >= Pstar)
    {
        result = Ucenter - (Ucenter & Normal) * Normal; // Тангенсальная компонента скорости
    }else
    {
        F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * MolWeightStar / Foam::constant::thermodynamic::RR / Tstar;
        U *= Foam::sqrt(Pstar/Rhostar); // Выход из нормировки на *
	Foam::vector N(Normal.x(), Normal.y(), Foam::sqrt(Normal.x()*Normal.x() + Normal.y()*Normal.y()) / 2.0);
        result = N / Foam::mag(N) * U; // Отнормируем вектор направления и умножим на скорость
    }
    return result;
}

Foam::scalar RhoinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[8];
    const Foam::scalar & Pstar = scalarParameters[9];
    const Foam::scalar & Tstar = scalarParameters[10];
    //------------------------------------------------------------------------
    Foam::scalar rho;
    Foam::scalar F, X, U;
    if(Pcenter >= Pstar)
    {
        rho = Rhocenter;
    }else
    {
        F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * MolWeightStar / Foam::constant::thermodynamic::RR / Tstar;
        rho = -F/U;
        rho *= Rhostar; // Выход из нормировки на *
    }
    return rho;
}

Foam::scalar TinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
    const Foam::scalar & SstarS = scalarParameters[8];
    const Foam::scalar & Pstar = scalarParameters[9];
    const Foam::scalar & Tstar = scalarParameters[10];
    //------------------------------------------------------------------------
    Foam::scalar Rho = RhoinSideWall(scalarParameters, time,  MolWeightStar, MolWeight, Pcenter, Tcenter, Rhocenter, Ucenter, c, gamma);
    Foam::scalar P = PinSideWall(scalarParameters, time, MolWeightStar, MolWeight, Pcenter, Tcenter, Rhocenter, Ucenter, c, gamma);
    return P * MolWeight / Rho / Foam::constant::thermodynamic::RR;
}


//================================================================================================================
// Выходящие параметры
template<class TScalar> TScalar SQRT(TScalar X);
template<class TScalar> TScalar POWER(TScalar X, TScalar Y); // X^Y

#include "../solver/breakdown/sln.H"
Foam::scalar Pout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
return Pcenter;
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & Backbressure = scalarParameters[1];
    //------------------------------------------------------------------------
    if(time < timeWall) return Pcenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    Foam::scalar cBackbressure = Foam::sqrt(gamma*Backbressure/Rhocenter);
    Foam::scalar ccenter = Foam::sqrt(gamma*Pcenter/Rhocenter);
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, cBackbressure, Pcenter, Rhocenter, Ucenter.z(), gamma, ccenter, P, Rho, U, Ustar, W);
    return P;
}

Foam::scalar Tout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
return Pcenter * MolWeight / Rhocenter / Foam::constant::thermodynamic::RR;
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & Backbressure = scalarParameters[1];
    //------------------------------------------------------------------------
    if(time < timeWall) return Pcenter * MolWeight / Rhocenter / Foam::constant::thermodynamic::RR;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    Foam::scalar cBackbressure = Foam::sqrt(gamma*Backbressure/Rhocenter);
    Foam::scalar ccenter = Foam::sqrt(gamma*Pcenter/Rhocenter);
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, cBackbressure, Pcenter, Rhocenter, Ucenter.z(), gamma, ccenter, P, Rho, U, Ustar, W);
    return P * MolWeight / Rho / Foam::constant::thermodynamic::RR;
}

Foam::vector Uout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
return Ucenter;
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & Backbressure = scalarParameters[1];
    //------------------------------------------------------------------------
    if(time < timeWall) return Ucenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    Foam::scalar cBackbressure = Foam::sqrt(gamma*Backbressure/Rhocenter);
    Foam::scalar ccenter = Foam::sqrt(gamma*Pcenter/Rhocenter);
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, cBackbressure, Pcenter, Rhocenter, Ucenter.z(), gamma, ccenter, P, Rho, U, Ustar, W);
    return Foam::vector(Ucenter.x(), Ucenter.y(), U);
}

Foam::scalar Rhoout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma)
{
return Rhocenter;
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & Backbressure = scalarParameters[1];
    //------------------------------------------------------------------------
    if(time < timeWall) return Rhocenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    Foam::scalar cBackbressure = Foam::sqrt(gamma*Backbressure/Rhocenter);
    Foam::scalar ccenter = Foam::sqrt(gamma*Pcenter/Rhocenter);
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, cBackbressure, Pcenter, Rhocenter, Ucenter.z(), gamma, ccenter, P, Rho, U, Ustar, W);
    return Rho;
}
