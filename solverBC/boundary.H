
#ifndef InletDictionary_H
#define InletDictionary_H

#include "scalarIOList.H"

Foam::scalar PinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar TinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter,  const Foam::scalar Tcenter,const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::vector UinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::vector Normal, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter,  const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar RhoinButt(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);

Foam::scalar PinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar TinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::vector UinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::vector Normal, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar RhoinSideWall(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);

Foam::scalar Pout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar Tout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::vector Uout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);
Foam::scalar Rhoout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar MolWeightStar, const Foam::scalar MolWeight, const Foam::scalar Pcenter, const Foam::scalar Tcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar c, const Foam::scalar gamma);

#endif