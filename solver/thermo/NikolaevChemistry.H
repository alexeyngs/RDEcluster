template<class ThermoType>
inline Foam::scalar GetMolWeightMaxNikolaev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Sigma)
{
    // Количество молей в смеси
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    // Количество атомов в смеси
    const Foam::scalar nH = 2.0*molesH2;
    const Foam::scalar nO = 2.0*molesO2;
    const Foam::scalar nN = 2.0*molesN2;
    const Foam::scalar nAr = molesAr;
    Foam::scalar MolWeightMax = 0.0;
    if(nO < 2.0*nH) MolWeightMax = Sigma / (nH/2.0 + nN/2.0 + nAr);
    if(nO >= nH/2.0) MolWeightMax = Sigma / (nH/4.0 + nO/2.0 + nN/2.0 + nAr);
    //--------------------------------------------------------------
    return MolWeightMax;
}

template<class ThermoType>
inline Foam::scalar GetQMaxNikolaev(const speciesTable & species, const PtrList<ThermoType> & speciesData)
{
    // Количество молей в смеси
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    // Вычисление энергии диссоциации Ed
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    Foam::scalar Qmax = 0.0;
    Foam::scalar Q0 = molesO2*Ed_O2 + molesH2*Ed_H2; // Отсчет от mu0
    if(molesO2 >= molesH2/2.0)   // избыток кислорода
    {
        Qmax = molesH2*Ed_H2O + (molesO2 - molesH2/2.0)*Ed_O2;
    }else // (molesO2 < molesH2/2.0) - недостаток кислорода
    {
        Qmax = 2.0*molesO2*Ed_H2O + (molesH2 - 2.0*molesO2)*Ed_H2;
    }
    Qmax = -Qmax;
	return Qmax;
}

template<class ThermoType>
inline Foam::scalar GetEdissociationNikolaev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Sigma, const Foam::scalar MolWeightMin)
{
    // Количество молей в смеси
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    // Вычисление энергии диссоциации Ed
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    Foam::scalar Qmax = 0.0;
    Foam::scalar Q0 = molesO2*Ed_O2 + molesH2*Ed_H2; // Отсчет от mu0
    if(molesO2 >= molesH2/2.0)   // избыток кислорода
    {
        Qmax = molesH2*Ed_H2O + (molesO2 - molesH2/2.0)*Ed_O2;
    }else // (molesO2 < molesH2/2.0) - недостаток кислорода
    {
        Qmax = 2.0*molesO2*Ed_H2O + (molesH2 - 2.0*molesO2)*Ed_H2;
    }
    Qmax = -Qmax;
    Q0 = -Q0;
    Foam::Info << "molesO2 = " << molesO2 << ", Qmax = " << Qmax/KiloCalorie << "Kcal (Qmax-Q0 = " << (Qmax-Q0)/KiloCalorie << "Kcal), Q0 = " << Q0/KiloCalorie << "Kcal" << Foam::endl;
    // Энергия диссоциации от min до max:
    Foam::scalar MolWeightMax = GetMolWeightMaxNikolaev(species, speciesData, Sigma);
    Foam::scalar Edissociation = Qmax / Sigma / (1.0/MolWeightMin - 1.0/MolWeightMax);
    return Edissociation;
}

template<class ThermoType>
inline Foam::scalar GetThetaNikolaev(const speciesTable & species, const PtrList<ThermoType> & speciesData)
{
    // Количество молей в смеси
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    // Энергия диссоциации Ed
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");

    const Foam::scalar theta_H2 = GetTheta(species, speciesData, "H2");
    const Foam::scalar theta_O2 = GetTheta(species, speciesData, "O2");
    const Foam::scalar theta_H2O = GetTheta(species, speciesData, "H2O");
    const Foam::scalar theta_N2 = GetTheta(species, speciesData, "N2");

    Foam::scalar theta = 0.0;
    if(molesO2 >= molesH2/2.0)   // избыток кислорода
    {
        theta = molesH2*theta_H2O + (molesO2 - molesH2/2.0)*theta_O2 + molesN2*theta_N2;
        theta /= (molesH2 + (molesO2 - molesH2/2.0) + molesN2);
    }else // (molesO2 < molesH2/2.0) - недостаток кислорода
    {
        theta = 2.0*molesO2*theta_H2O + (molesH2 - 2.0*molesO2)*theta_H2 + molesN2*theta_N2;
        theta /= 2.0*molesO2 + (molesH2 - 2.0*molesO2) + molesN2;   
    }
    return theta;
}