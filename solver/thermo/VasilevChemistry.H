template<class ThermoType>
inline Foam::scalar GetMolWeightFromInductionVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Induction, const Foam::scalar powerN, const Foam::scalar Sigma)
{
    const Foam::scalar BETA = Foam::pow(Induction, powerN);
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    //--------------------------------------------------------------
    Foam::scalar MolWeight = Sigma / (2.5*molesCH4 - 1.5*BETA*molesCH4 + molesH2 + molesO2 + molesN2 + molesAr + molesH2O);
    return MolWeight;
}

template<class ThermoType>
inline Foam::scalar GetInductionFromHchemicalVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Hchemical, const Foam::scalar powerN, const Foam::scalar Sigma)
{
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = GetEdissociation(species, speciesData, "CH4");
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    const Foam::scalar Ed_OH = GetEdissociation(species, speciesData, "OH");
    //--------------------------------------------------------------
    Foam::scalar H = Hchemical / 1000.0 * Sigma;
    Foam::scalar NoBeta = molesCH4*Ed_OH + 1.5*molesCH4*Ed_H2 + molesH2*Ed_H2 + molesO2*Ed_O2 - molesCH4*Ed_O2 + molesH2O*Ed_H2O;
    Foam::scalar WithBeta = molesCH4*Ed_CH4 - molesCH4*Ed_OH - 1.5*molesCH4*Ed_H2 + molesCH4*Ed_O2;
    Foam::scalar BETA = (H - NoBeta) / WithBeta;
    if(BETA > 1.0-1.0e-12) BETA = 1.0;
    if(BETA < 0.0) return 0.0;
    Foam::scalar Induction = Foam::pow(BETA, 1.0/powerN);
    if(Induction > 1.0) return Induction = 1.0;
    return Induction;
}

template<class ThermoType>
inline Foam::scalar GetHchemicalFromInductionVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Induction, const Foam::scalar powerN, const Foam::scalar Sigma)
{
    const Foam::scalar BETA = Foam::pow(Induction, powerN);
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = GetEdissociation(species, speciesData, "CH4");
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    const Foam::scalar Ed_OH = GetEdissociation(species, speciesData, "OH");
    Foam::scalar Hchemical = (BETA*molesCH4*Ed_CH4 + (1.0 - BETA)*molesCH4*Ed_OH + (1.5*molesCH4 - 1.5*BETA*molesCH4 + molesH2)*Ed_H2 + (molesO2 - molesCH4 + BETA*molesCH4)*Ed_O2 + molesH2O*Ed_H2O) / Sigma;
    return Hchemical*1000.0;
}

template<class ThermoType>
inline Foam::scalar GetMolWeightMaxVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Sigma)
{
    // Количество молей в смеси
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    // Количество атомов в смеси
    const Foam::scalar nC = molesCH4;
    const Foam::scalar nH = 4.0*molesCH4 + 2.0*molesH2 + 2.0*molesH2O;
    const Foam::scalar nO = 2.0*molesO2 + molesH2O;
    const Foam::scalar nN = 2.0*molesN2;
    const Foam::scalar nAr = molesAr;
    //--------------------------------------------------------------
    Foam::scalar MolWeightMax = 0.0;
    if(nC <= nO < 2.0*nC) MolWeightMax = Sigma / (nC + nH/2.0 + nN/2.0 + nAr);
    if(2.0*nC <= nO < 2.0*nH) MolWeightMax = Sigma / (nC + nH/2.0 + nN/2.0 + nAr);
    if(nO >= 2.0*nC + nH/2.0) MolWeightMax = Sigma / (nH/4.0 + nO/2.0 + nN/2.0 + nAr);
    return MolWeightMax;
}

template<class ThermoType>
inline Foam::scalar GetQMaxVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData)
{
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = GetEdissociation(species, speciesData, "CH4");
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    const Foam::scalar Ed_CO2 = GetEdissociation(species, speciesData, "CO2");
    const Foam::scalar Ed_OH = GetEdissociation(species, speciesData, "OH");
    const Foam::scalar Ed_CO = 0.0; // Точка отсчета от CO, поэтому EdCO равен 0
    Foam::scalar Qmax = 0.0;
    Foam::scalar Qstar = 0.0;
    //--------------------------------------------------------------
    if(molesCH4/2.0 <= molesO2 && molesO2 < molesCH4)
    {
        // Кислорода хватило только на окисление углерода и немного водорода
        // Разложение топлива CH4 + O2 -> CO + H2
        Qstar = (2.0*molesO2 - molesCH4)*Ed_CO + (2.0*molesO2 - molesCH4)*Ed_OH + (2.5*molesCH4 + molesH2 - molesO2)*Ed_H2;
    }else
    {
        // Кислорода хватит полностью на CO и ОН
        // Разложение топлива CH4 + O2 -> CO + OH + 1.5H2
        Qstar = molesCH4*Ed_CO + molesCH4*Ed_OH + 1.5*molesCH4*Ed_H2 + (molesO2 - molesCH4)*Ed_O2 + molesH2*Ed_H2;
    }
    //--------------------------------------------------------------
    if(molesO2 >= 2.0*molesCH4 + molesH2/2.0)   // Явно избыток кислорода
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesCH4)*Ed_H2O + molesH2*Ed_H2O + (molesO2 - 2.0*molesCH4 - molesH2/2.0)*Ed_O2;
    }else if(molesCH4 <= molesO2 && molesO2 < 2.0*molesCH4 + molesH2/2.0)   // Кислорода хватит только на окисления углерода CO2
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesO2 - 2.0*molesCH4)*Ed_H2O + (4.0*molesCH4 + molesH2 - 2.0*molesO2)*Ed_H2;
    }else if(molesO2 < molesCH4)    // Кислорода не хватило даже на полное окисление углерода, но на образование CO кислорода есть всегда, т.к. должны получать газ
    {
        Qmax = (2.0*molesO2 - molesCH4)*Ed_CO2 + 2.0*molesCH4*Ed_H2;
    };
    Foam::scalar Q0 = molesCH4*Ed_CH4 + molesO2*Ed_O2 + molesH2*Ed_H2; // Отсчет от mu0
    Qmax = -Qmax;
    return Qmax;
}

template<class ThermoType>
inline Foam::scalar GetEdissociationVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Sigma, const Foam::scalar MolWeightMin)
{
    const Foam::scalar molesCH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar molesH2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar molesO2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar molesN2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar molesAr = GetMoles(species, speciesData, "Ar");
    const Foam::scalar molesH2O = GetMoles(species, speciesData, "H2O");
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = GetEdissociation(species, speciesData, "CH4");
    const Foam::scalar Ed_H2 = GetEdissociation(species, speciesData, "H2");
    const Foam::scalar Ed_O2 = GetEdissociation(species, speciesData, "O2");
    const Foam::scalar Ed_H2O = GetEdissociation(species, speciesData, "H2O");
    const Foam::scalar Ed_CO2 = GetEdissociation(species, speciesData, "CO2");
    const Foam::scalar Ed_OH = GetEdissociation(species, speciesData, "OH");
    const Foam::scalar Ed_CO = 0.0; // Точка отсчета от CO, поэтому EdCO равен 0
    Foam::scalar Qmax = 0.0;
    Foam::scalar Qstar = 0.0;
    //--------------------------------------------------------------
    if(molesCH4/2.0 <= molesO2 && molesO2 < molesCH4)
    {
        // Кислорода хватило только на окисление углерода и немного водорода
        // Разложение топлива CH4 + O2 -> CO + H2
        Qstar = (2.0*molesO2 - molesCH4)*Ed_CO + (2.0*molesO2 - molesCH4)*Ed_OH + (2.5*molesCH4 + molesH2 - molesO2)*Ed_H2;
    }else
    {
        // Кислорода хватит полностью на CO и ОН
        // Разложение топлива CH4 + O2 -> CO + OH + 1.5H2
        Qstar = molesCH4*Ed_CO + molesCH4*Ed_OH + 1.5*molesCH4*Ed_H2 + (molesO2 - molesCH4)*Ed_O2 + molesH2*Ed_H2;
    }
    //--------------------------------------------------------------
    if(molesO2 >= 2.0*molesCH4 + molesH2/2.0)   // Явно избыток кислорода
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesCH4)*Ed_H2O + molesH2*Ed_H2O + (molesO2 - 2.0*molesCH4 - molesH2/2.0)*Ed_O2;
    }else if(molesCH4 <= molesO2 && molesO2 < 2.0*molesCH4 + molesH2/2.0)   // Кислорода хватит только на окисления углерода CO2
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesO2 - 2.0*molesCH4)*Ed_H2O + (4.0*molesCH4 + molesH2 - 2.0*molesO2)*Ed_H2;
    }else if(molesO2 < molesCH4)    // Кислорода не хватило даже на полное окисление углерода, но на образование CO кислорода есть всегда, т.к. должны получать газ
    {
        Qmax = (2.0*molesO2 - molesCH4)*Ed_CO2 + 2.0*molesCH4*Ed_H2;
    };
    Foam::scalar Q0 = molesCH4*Ed_CH4 + molesO2*Ed_O2 + molesH2*Ed_H2; // Отсчет от mu0
    Qmax = -Qmax;
    Q0 = -Q0;
    Qstar = -Qstar;
    Info << "molesO2 = " << molesO2 << ", Qmax = " << Qmax/4.1868 << "(" << (Qmax-Q0)/4.1868<< ")" << " cal, Q* = " << Qstar/4.1868 << "(" << (Qstar-Q0)/4.1868<< ")" << "cal , Q0 = " << Q0/4.1868 << " cal" << endl;
    Foam::scalar MolWeightStar = GetMolWeightFromInduction(species, speciesData, 0.0, 1.0, Sigma);
    Foam::scalar MolWeightMax = GetMolWeightMaxVasilev(species, speciesData, Sigma);
    // Энергия диссоциации от * до max:
//    Foam::scalar Edissociation = (Qmax - Qstar) / Sigma / (1.0/MolWeightStar - 1.0/MolWeightMax);
    // Энергия диссоциации от min до max:
    Foam::scalar Edissociation = Qmax / Sigma / (1.0/MolWeightMin - 1.0/MolWeightMax);
    return Edissociation;
}

template<class ThermoType>
inline Foam::scalar GetThetaVasilev(const speciesTable & species, const PtrList<ThermoType> & speciesData)
{
    const Foam::scalar moles_CH4 = GetMoles(species, speciesData, "CH4");
    const Foam::scalar moles_H2 = GetMoles(species, speciesData, "H2");
    const Foam::scalar moles_O2 = GetMoles(species, speciesData, "O2");
    const Foam::scalar moles_N2 = GetMoles(species, speciesData, "N2");
    const Foam::scalar moles_H2O = GetMoles(species, speciesData, "H2O");
    const Foam::scalar moles_CO2 = GetMoles(species, speciesData, "CO2");
    const Foam::scalar moles_OH = GetMoles(species, speciesData, "OH");
    const Foam::scalar moles_CO = GetMoles(species, speciesData, "OH");
    //--------------------------------------------------------------
    const Foam::scalar theta_CH4 = GetTheta(species, speciesData, "CH4");
    const Foam::scalar theta_H2 = GetTheta(species, speciesData, "H2");
    const Foam::scalar theta_O2 = GetTheta(species, speciesData, "O2");
    const Foam::scalar theta_N2 = GetTheta(species, speciesData, "N2");
    const Foam::scalar theta_H2O = GetTheta(species, speciesData, "H2O");
    const Foam::scalar theta_CO2 = GetTheta(species, speciesData, "CO2");
    const Foam::scalar theta_OH = GetTheta(species, speciesData, "OH");
    const Foam::scalar theta_CO = GetTheta(species, speciesData, "OH");
    //--------------------------------------------------------------
    Foam::scalar theta = 0.0;
    // if(moles_CH4/2.0 <= moles_O2 && moles_O2 < moles_CH4)
    // {
    //     theta = (2.0*moles_O2 - moles_CH4)*theta_CO + (2.0*moles_O2 - moles_CH4)*theta_OH + (2.5*moles_CH4 + moles_H2 - moles_O2)*theta_H2 + moles_N2*theta_N2;
    //     theta /= (2.0*moles_O2 - moles_CH4) + (2.0*moles_O2 - moles_CH4) + (2.5*moles_CH4 + moles_H2 - moles_O2) + moles_N2;
    // }else
    // {
    //     theta = moles_CH4*theta_CO + moles_CH4*theta_OH + 1.5*moles_CH4*theta_H2 + (moles_O2 - moles_CH4)*theta_O2 + moles_H2*theta_H2 + moles_N2*theta_N2;
    //     theta /= moles_CH4 + moles_CH4 + 1.5*moles_CH4 + (moles_O2 - moles_CH4) + moles_H2 + moles_N2;
    // }
    if(moles_O2 >= 2.0*moles_CH4 + moles_H2/2.0)   // Явно избыток кислорода
    {
        theta = moles_CH4*theta_CO2 + (2.0*moles_CH4)*theta_H2O + moles_H2*theta_H2O + (moles_O2 - 2.0*moles_CH4 - moles_H2/2.0)*theta_O2 + moles_N2*theta_N2;
        theta /= moles_CH4 + (2.0*moles_CH4) + moles_H2 + (moles_O2 - 2.0*moles_CH4 - moles_H2/2.0) + moles_N2;
    }else if(moles_CH4 <= moles_O2 && moles_O2 < 2.0*moles_CH4 + moles_H2/2.0)   // Кислорода хватит только на окисления углерода CO2
    {
        theta = moles_CH4*theta_CO2 + (2.0*moles_O2 - 2.0*moles_CH4)*theta_H2O + (4.0*moles_CH4 + moles_H2 - 2.0*moles_O2)*theta_H2 + moles_N2*theta_N2;
        theta /= moles_CH4 + (2.0*moles_O2 - 2.0*moles_CH4) + (4.0*moles_CH4 + moles_H2 - 2.0*moles_O2) + moles_N2;
    }else if(moles_O2 < moles_CH4)    // Кислорода не хватило даже на полное окисление углерода, но на образование CO кислорода есть всегда, т.к. должны получать газ
    {
        theta = (2.0*moles_O2 - moles_CH4)*theta_CO2 + 2.0*moles_CH4*theta_H2 + moles_N2*theta_N2;
        theta /= (2.0*moles_O2 - moles_CH4) + 2.0*moles_CH4 + moles_N2;
    };
    return theta;
}