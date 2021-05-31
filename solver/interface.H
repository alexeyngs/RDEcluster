#ifndef INTERFACE_H
#define INTERFACE_H
namespace Foam
{
    struct IMixture
    {
        virtual scalar GetA(const scalar MolWeight, const scalar T) const = 0;
        virtual scalar Cv(const scalar P, const scalar T, const scalar MolWeight) const = 0;
        virtual scalar GetTmax(const scalar MolWeight, const scalar Uthermodinamical) const = 0;
        virtual const volScalarField & Rho() const = 0;
        virtual const volScalarField & he() const = 0;
    };
};

enum class EBreakDown
{
    None,
    Accustic,
    Godunov,
    Kolgan,
};

enum class EUpLayer
{
    THE,
    getT,
};

const Foam::scalar Epsilon = Foam::SMALL*10.0;
#endif //INTERFACE_H