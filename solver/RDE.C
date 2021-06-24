const unsigned int NVersion = 85;
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "pointPatch.H"
#include "scalarIOList.H"

#include "interface.H"
#include "thermo/rhoRDEThermo.H"
//Foam::constant::thermodynamic

Foam::word InfoBreakDown();
void Stitch(Foam::Time & runTime, Foam::polyMesh & mesh, const Foam::word masterPatchName, const Foam::word slavePatchName);
void ReadBreakDown(const Foam::dictionary & dict);
void BreakDown(Foam::rhoRDEThermo & Thermo, const Foam::volScalarField & Psi, const Foam::volScalarField & c, const Foam::volScalarField & gamma, const Foam::volScalarField & T,
    const Foam::volScalarField & Rho, const Foam::volVectorField & U, const Foam::volScalarField & P, const Foam::volScalarField & E, const Foam::volScalarField & MolWeight, const Foam::volScalarField & Induction,
    const Foam::surfaceScalarField & Rholf, const Foam::surfaceVectorField & Ulf, const Foam::surfaceScalarField & Plf, const Foam::surfaceScalarField & Elf, const Foam::surfaceScalarField & MolWeightlf, const Foam::surfaceScalarField & Inductionlf,
    Foam::surfaceScalarField & Rhof, Foam::surfaceVectorField & Uf, Foam::surfaceScalarField & Pf, Foam::surfaceScalarField & Ef, Foam::surfaceScalarField & Hthermodynamicalf, Foam::surfaceScalarField & MolWeightf, Foam::surfaceScalarField & Inductionf);

// Работа с жесткой стенкой
#include "solidWall.H"
//----------------------------------------------------------------------------------------------------------------
bool TimeWall = true;
bool EmptyContactMesh = false;
Foam::label Multiplier = 1;
Foam::scalar CameraLength = 0.0;
Foam::scalar CameraDiameter = 0.0;
Foam::scalar Tmin = 0.0;
EUpLayer UpLayer = EUpLayer::THE;
int main(int argc, char * argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createParameters.H"
    #include "createFields.H"
    #include "createFace.H"
    #include "createTimeControls.H"
    //----------------------------------------------------------------------------------------------------------------
    #include "contact/createContactFields.H"
    // Courant numbers used to adjust the time-step
    Courant = (mag(U) + c)*runTime.deltaT() / hdelta;
    Foam::scalar CoNum = max(Courant).value();
    // Инициализация переключателя на периодические условия
    timeOnePerimeter = CameraDiameter * Foam::constant::mathematical::pi / thermo.Dcj();
    TimeWall = runTime.value() < timeOnePerimeter * TimeKoefficientDetonationPerimeter;
    bool LastTimeSolidWall = TimeWall;
    //----------------------------------------------------------------------------------------------------------------
    Foam::Info << "Count processors " << Foam::Pstream::nProcs() << Foam::endl;
    Foam::Info << "My rank in MPI Communicator is " << Pstream::myProcNo() << " and master rank " << Pstream::masterNo() << Foam::endl;
    Foam::Info << "timeOnePerimeter = " << timeOnePerimeter << Foam::endl;
    Foam::label iter = runTime.startTimeIndex();
    Foam::Info << "Start iterator = " << iter << Foam::endl;
    Foam::Info << "Count cells = " << MolWeight.size() << Foam::endl;
    Foam::Info << "Size long double = " << sizeof(long double) << Foam::endl;
    Foam::Info << "CameraLength = " << CameraLength << Foam::endl;
    Foam::Info << "CameraDiameter = " << CameraDiameter << Foam::endl;
    Foam::Info << "UseChemistry = " << UseChemistry << Foam::endl;
    Foam::Info << Backpressure << Foam::endl;
    Foam::Info << "TimeKoefficientDetonationPerimeter = " << TimeKoefficientDetonationPerimeter << Foam::endl;
    Foam::Info << "Tmin = " << Tmin << "K" << Foam::endl;
    Foam::Info << "deltaTmin = " << deltaTmin << " second" << Foam::endl;
    Foam::Info << "LastTimeSolidWall = " << LastTimeSolidWall << Foam::endl;
    Foam::Info << "Multiplier = " << Multiplier << Foam::endl;
    Foam::Info << "BreakDown: " << InfoBreakDown() << Foam::endl;
    Foam::Info << "Tuplayer: ";
    switch(UpLayer)
    {
        case EUpLayer::THE : Foam::Info << "THE"; break;
        case EUpLayer::getT : Foam::Info << "getT"; break;
        default: Foam::Info << "error!"; break;
    }Foam::Info << Foam::endl;
    Foam::Info << "RDE version-" << NVersion << Foam::endl;
    Foam::Info << "Starting time LOOP!" << Foam::endl;
    Foam::Info << "==========================================================================================================================" << endl << endl;
    if(UseChemistry) thermo.chemistry(CameraLength); // Один раз вызовем химию для заполнения правой части
    while (runTime.run())
    {
        BreakDown(thermo, psi, c, gamma, T, rho, U, p, e, MolWeight, Induction, Rholf, Ulf, Plf, Elf, MolWeightlf, Inductionlf, Rhof, Uf, Pf, Ef, Hthermodynamicalf, MolWeightf, Inductionf);
        thermo.Paverage = fvc::average(Pf);
        thermo.Rhoaverage = fvc::average(Rhof);
        thermo.Paverage.correctBoundaryConditions();
        thermo.Rhoaverage.correctBoundaryConditions();
        // Временная жесткая стенка
        if(TimeWall)
        {
            SolidWall(U, Uf, p, Pf, MolWeight, MolWeightf, Induction, Inductionf);
        }else if(LastTimeSolidWall)
        {
            LastTimeSolidWall = false;
            LastSolidWall(e, T, U, p, MolWeight, Induction);
            U.correctBoundaryConditions();
            e.correctBoundaryConditions();
            p.correctBoundaryConditions();
            MolWeight.correctBoundaryConditions();
            Induction.correctBoundaryConditions();
            rho = thermo.CorrectRho();
            rhoU = rho*U;
            rhoE = rho*(e + 0.5*magSqr(U));
            rhoInduction = rho*Induction;
            rhoMolWeight = rho*MolWeight;
            // ребра при линейной интерполяции
            Rholf = linearInterpolate(rho);
            Ulf = linearInterpolate(U);
            Plf = linearInterpolate(p);
            Tlf = linearInterpolate(T);
            Elf = linearInterpolate(e);
            MolWeightlf = linearInterpolate(MolWeight);
            Inductionlf = linearInterpolate(Induction);
            // Нахождение температуры
            thermo.correct();
            // Исправление ошибок по температуре
            thermo.CorrectErrors(UseChemistry, rho, Rhof, Pf);
            gamma.ref() = thermo.Cp()/thermo.Cv();
            // Скорость звука
            c.ref() = Foam::sqrt(gamma / psi);
            BreakDown(thermo, psi, c, gamma, T, rho, U, p, e, MolWeight, Induction, Rholf, Ulf, Plf, Elf, MolWeightlf, Inductionlf, Rhof, Uf, Pf, Ef, Hthermodynamicalf, MolWeightf, Inductionf);
            thermo.Paverage = fvc::average(Pf);
            thermo.Rhoaverage = fvc::average(Rhof);
            thermo.Paverage.correctBoundaryConditions();
            thermo.Rhoaverage.correctBoundaryConditions();
        }

        #include "updateFace.H"
        #include "readTimeControls.H"
        //<<<<<<<<<<<<<<<<<<<<<<<<
        Foam::Info << nl << "-------------------------" << Foam::endl;
        #include "setDeltaT.H"
        if(runTime.deltaT().value() < deltaTmin)
        {
            runTime.setDeltaT(deltaTmin);
            Foam::Info << "deltaT corrected, = " <<  runTime.deltaTValue() << Foam::endl;
        }
        Courant = (mag(U) + c)*runTime.deltaT() / hdelta;
        Foam::Info << "begin" << Foam::endl;
        CoNum = max(Courant).value();
        Foam::Info << "end" << Foam::endl;
        Foam::Info << "Courant Number = " << CoNum << Foam::endl;
        Foam::Info << "Iterator = " << iter << Foam::endl;
        //>>>>>>>>>>>>>>>>>>>>>>>>
        // Переключение от жесткой стенки к периодическим условиям
        TimeWall = runTime.value() < timeOnePerimeter * TimeKoefficientDetonationPerimeter;
        runTime++;
        iter++;
        Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

        // Уравнение массы
        solve
        (
            fvm::ddt(rho) + fvc::div(RhoUf)
        );

        // Уравнение импульса
        solve
        (
            fvm::ddt(rhoU) + fvc::div(RhoUUf) + fvc::div(PPf)
        );

        // Вычисление скорости
        U.ref() = rhoU() / rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();


        // Уравнение энергии
        solve
        (
            fvm::ddt(rhoE) + fvc::div(RhoUHf)
        );
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5*magSqr(U.boundaryField()));

	    // Химия
        if(UseChemistry) thermo.chemistry(CameraLength);
        
        // Уравнение доли периода индукции Induction
        solve
        (
            fvm::ddt(rhoInduction) + fvc::div(RhoUInductionf) - FInduction
        );
        Induction.ref() = rhoInduction() / rho();
        forAll(Induction, i)
        {
            if(Induction[i] < 0.0) Induction[i] = 0.0;
            //if(Foam::mag(Zmax[i] - CameraLength) < Epsilon) Induction[i] = InductionZmax[i];
        }
        Induction.correctBoundaryConditions();
        rhoInduction.ref() = rho*Induction;
        rhoInduction.boundaryFieldRef() == rho.boundaryField()*Induction.boundaryField();

        // Уравнение молярной массы MolWeight
        solve
        (
            fvm::ddt(rhoMolWeight) + fvc::div(RhoUMolWeightf) - FMolWeight
        );
        // Подготовка молярной массы
        MolWeight.ref() = rhoMolWeight() / rho();
        thermo.CorrectChemistry();
        // Контактная граница
        #include "contact/updateContactFields.H"
        MolWeight.correctBoundaryConditions();
thermo.CorrectChemistry();
        // Конец создания молярной массы
        rhoMolWeight.ref() = rho*MolWeight;
        rhoMolWeight.boundaryFieldRef() == rho.boundaryField()*MolWeight.boundaryField();
        //----------------------------------------------------------------------
        // Нахождение температуры
        thermo.correct();
        // Исправление ошибок по температуре
        thermo.CorrectErrors(UseChemistry, rho, Rhof, Pf);
        gamma.ref() = thermo.Cp()/thermo.Cv();
        // Скорость звука
        c.ref() = sqrt(gamma / psi);
        // Число Маха вдоль оси Z
        Mz = -U.component(vector::Z) / c;

        // Correct pressure
        p.ref() = rho() / psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        g = -rho*U.component(vector::Z);
        impuls = p + rho*U.component(vector::Z)*U.component(vector::Z) - Backpressure;
        forAll(Ucylinder, i)
        {
            const Foam::vector Centre = mesh.C()[i]; // Координата центра ячейки
            const Foam::scalar Ur = (U[i].x()*Centre.x() + U[i].y()*Centre.y()) / Foam::sqrt(Centre.x()*Centre.x() + Centre.y()*Centre.y());
            // phi. = (xy. - x.y)/(xx+yy)
            // Uphi = r*phi.
            Foam::scalar Uphi = (Centre.x()*U[i].y() - Centre.y()*U[i].x()) / Foam::sqrt(Centre.x()*Centre.x() + Centre.y()*Centre.y());
            Ucylinder[i] = vector(Ur, Uphi, U[i].z());
        }
        runTime.write();
        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << Foam::endl;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//#include "file/FieldFile.H"
//#include "file/ContactFile.H"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        forAll(p, celli)
        {
            if(T[celli] < 0.0 || p[celli] < 0.0 || rho[celli] < 0.0)
            {
                Foam::scalar PP = p[celli];
                Foam::scalar TT = T[celli];
                Foam::scalar rhorho = rho[celli];
                Foam::Info << "CONTROL ERROR! celli = " << celli << " P = " << PP << " , T = " << TT << " , rho = " << rhorho << Foam::endl;
            }
        }
    }
    Foam::Info << "End" << Foam::endl;
    return 0;
}

template<typename TScalar>
TScalar SQRT(TScalar X);

template<typename TScalar>
TScalar POWER(TScalar X, TScalar Y); // X^Y
//====================================================================================================================
template<>
Foam::scalar SQRT<Foam::scalar>(Foam::scalar X)
{
    return Foam::sqrt(X);
};

template<>
Foam::scalar POWER<Foam::scalar>(Foam::scalar X, Foam::scalar Y) // X^Y
{
	return Foam::pow(X,Y);
};

// Основные положения:
// Главная ось симметрии - обязательно ось Z.
// Камера сгорания имеет диапазон по оси Z от величины CameraLength(вход) до 0(выход).
// Всё, что выше по оси Z - CameraLength - это область вне камеры сгорания, поэтому химия там не включается.
// Химия работает только в ячейках, у которых значение Z ниже значения параметра CameraLength.
// Параметры CameraLength и CameraDiameter устанавливается в файле controlDict.


// ХИМИЧЕСКАЯ КИНЕТИКА
// Николаев - для водородо воздушной смеси (H2 + air) и для синтез-газа (CO + H2 + air)
// Васильев - для метано-водородо воздушной смеси (CH4 + H2 + air)

