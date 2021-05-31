#include "rhoRDEThermo.H"
#include "fvMesh.H"
#include "specie.H"
#include "IOobject.H"


#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
namespace Foam
{
    defineTypeNameAndDebug(rhoRDEThermo, 0);
    defineRunTimeSelectionTable(rhoRDEThermo, fvMesh);
}
Foam::rhoRDEThermo::rhoRDEThermo(const fvMesh& mesh,const word& phaseName) : Foam::rhoReactionThermo(mesh, phaseName),
Paverage
(
    IOobject
    (
        "Paverage",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::average(linearInterpolate(this->p_))
),
Rhoaverage
(
    IOobject
    (
        "Rhoaverage",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::average(linearInterpolate(this->rho_))
)
{
}
Foam::autoPtr<Foam::rhoRDEThermo> Foam::rhoRDEThermo::New(const fvMesh& mesh,const word& phaseName)
{
    return basicThermo::New<rhoRDEThermo>(mesh, phaseName);
}
Foam::rhoRDEThermo::~rhoRDEThermo()
{
}
