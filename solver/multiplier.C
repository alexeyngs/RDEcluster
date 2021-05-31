#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Field.H"
#include "thermo/rhoRDEThermo.H"
#include <vector>
#include <utility>


typedef std::pair<Foam::scalar, Foam::Field<Foam::label> > PANCAKE;
Foam::Field<Foam::label> MultiplyBase;
inline Foam::scalar Angle(const Foam::vector & Position)
{
    Foam::scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
    return Position.y() > 0 ? phi : 2.0*Foam::constant::mathematical::pi - phi;
};

void CreateMultiplyBase(const Foam::label Multiplier, const Foam::fvMesh & mesh, const Foam::scalar CameraLength, const Foam::scalar MinHdelta)
{
    const Foam::scalar Epsilon = Foam::SMALL*100.0;
    // Создаем размер
    MultiplyBase.resize(mesh.nCells());
    // Список блинов
    std::vector<PANCAKE> Pancakes;
    forAll(mesh.C(), index)
    {
        const Foam::vector Position = mesh.C()[index];
        MultiplyBase[index] = index;
        // Если индекс находится в камере сгорания
        if(Position.z() < CameraLength + Epsilon)
        {
            // Сортируем в блины
            bool added = false;
            for(PANCAKE & elem : Pancakes) if(Foam::mag(Position.z() - elem.first) < Epsilon)
            {
                // Увеличиваем размер блина при определенном значения Z
                elem.second.append(index);
                added = true;
                break;
            }
            // Создаем новый блин
            if(!added)
            {
                PANCAKE cell;
                cell.first = Position.z();
                cell.second.append(index);
                Pancakes.push_back(cell);
            }
        }
    }

    // Разбираемся с блинами
    for(PANCAKE & elem : Pancakes)
    {
        forAll(elem.second, index)
        {
            const Foam::vector POSITION = mesh.C()[elem.second[index]];
            const Foam::scalar R = Foam::sqrt(POSITION.x() * POSITION.x() + POSITION.y() * POSITION.y());
            Foam::scalar PHI = Angle(POSITION);
            while(PHI > 2.0*Foam::constant::mathematical::pi/Multiplier) PHI -= 2.0*Foam::constant::mathematical::pi/Multiplier;
            PHI *= Multiplier;
            Foam::scalar delta = Foam::GREAT;
            // Поиск по углу
            forAll(elem.second, findindex)
            {
                Foam::vector position = mesh.C()[elem.second[findindex]];
                Foam::scalar r = Foam::sqrt(position.x() * position.x() + position.y() * position.y());
                if(Foam::mag(R - r) > MinHdelta) continue;
                if(Foam::mag(PHI - Angle(position)) < delta)
                {
                    delta = Foam::mag(PHI - Angle(position));
                    MultiplyBase[elem.second[index]] = elem.second[findindex];
                }
            }
            // Конец поиска по углу
        }
    }
}

void CreateMultiplyWave(Foam::volScalarField & Field)
{
    Foam::volScalarField TempField(Field);
    forAll(Field, i) Field[i] = TempField[MultiplyBase[i]];
    TempField.clear();
}
void CreateMultiplyWave(Foam::volVectorField & Field, const Foam::fvMesh & mesh)
{
    Foam::volVectorField TempField(Field);
    forAll(Field, i) if(i == MultiplyBase[i]) Field[i] = TempField[MultiplyBase[i]];
    else
    {
        // Считаем, что инвариантом будет угловая, радиальная, Z составлющая вектора(скорости)
        // То есть Fr, Fphi, Fz остаются неизменными при любом изменении угла phi
        Foam::vector Pold = mesh.C()[MultiplyBase[i]];  // Position - координата центра ячейки
        Foam::vector Pnew = mesh.C()[i];                // Position - координата центра ячейки
        Foam::vector F = TempField[MultiplyBase[i]];    // Field - ячейка векторного поля
        // Нораль:      (x, y)
        // Тангенсаль:  (-y,x)
        const Foam::scalar Fr = (Pold.x()*F.x() + Pold.y()*F.y()) / Foam::sqrt(Pold.x()*Pold.x() + Pold.y()*Pold.y());
        const Foam::scalar Ft = (-Pold.y()*F.x() + Pold.x()*F.y()) / Foam::sqrt(Pold.x()*Pold.x() + Pold.y()*Pold.y());
        // Заолнение нового вектора в ячейку
        Field[i] = Foam::vector
        (
            (Pnew.x()*Fr - Pnew.y()*Ft) / Foam::sqrt(Pnew.x()*Pnew.x() + Pnew.y()*Pnew.y()),
            (Pnew.y()*Fr + Pnew.x()*Ft) / Foam::sqrt(Pnew.x()*Pnew.x() + Pnew.y()*Pnew.y()),
            F.z()
        );
    }
    TempField.clear();
}

void MultiplyWave(const Foam::label Multiplier, const Foam::scalar CameraLength, const Foam::scalar MinHdelta, const Foam::fvMesh & mesh, Foam::rhoRDEThermo & thermo, Foam::volScalarField & rho, Foam::volVectorField & U)
{
    CreateMultiplyBase(Multiplier, mesh, CameraLength, MinHdelta);
    CreateMultiplyWave(rho);
    CreateMultiplyWave(U, mesh);
    CreateMultiplyWave(thermo.rho());
    CreateMultiplyWave(thermo.he());
    CreateMultiplyWave(thermo.p());
    CreateMultiplyWave(thermo.T());
    CreateMultiplyWave(const_cast<Foam::volScalarField &> (thermo.psi()));
    CreateMultiplyWave(thermo.MolWeight());
    CreateMultiplyWave(thermo.Induction());
}




