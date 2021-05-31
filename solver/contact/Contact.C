#include "IOmanip.H"
#include "fvCFD.H"
#include "pTraits.H"
#include "Field.H"
#include "boundaryMesh.H"
#include "thermo/rhoRDEThermo.H"
#include "../interface.H"
inline Foam::scalar Angle(const Foam::vector & Position)
{
    Foam::scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
    return Position.y() > 0 ? phi : 2.0*Foam::constant::mathematical::pi - phi;
};

// Функция получения номера ячейки на блине контактной границе в зависимости от номера ячейки в камере сгорания
inline Foam::label GetLabelBounaryZ(const Foam::fvMesh & mesh, const Foam::vectorField & ContactCoordinateMesh, const Foam::label celli)
{
    // Точка (X,Y) для конкретной ячейки
    const Foam::point Position(mesh.cellCentres()[celli].x(), mesh.cellCentres()[celli].y(), 0.0);
    Foam::label result = -1;
    Foam::scalar length = Foam::GREAT;
    forAll(ContactCoordinateMesh, i) // По всем ячейкам
    {
        Foam::point P(ContactCoordinateMesh[i].x(), ContactCoordinateMesh[i].y(), 0.0);
        Foam::scalar delta = Foam::mag(P - Position);
        if(length > delta)
        {
            length = delta;
            result = i;
        }
    }
    return result;
}

// Функция получения номера ячейки на поверхности контактной границе в зависимости от номера ячейки в камере сгорания
inline Foam::label GetLabelBounaryR(const Foam::fvMesh & mesh, const Foam::vectorField & ContactCoordinateMesh, const Foam::label celli)
{
    const Foam::scalar Epsilon = Foam::SMALL*10.0;
    // Точка для угла
    const Foam::point & POSITION = mesh.cellCentres()[celli];
    const Foam::scalar PHI = Angle(POSITION);
    Foam::label result = -1;
    Foam::scalar length = Foam::GREAT;
    forAll(ContactCoordinateMesh, i) // По всем ячейкам
    {
        const Foam::point & position = ContactCoordinateMesh[i];
        const Foam::scalar phi = Angle(position);
        Foam::scalar deltaPhi = Foam::mag(PHI - phi);
        Foam::scalar deltaZ = Foam::mag(POSITION.z() - position.z());
        if(length > deltaPhi && deltaZ < Epsilon)
        {
            length = deltaPhi;
            result = i;
        }
    }
    return result;
}

void UpdateContactMesh(
    Foam::volScalarField & Zmin, Foam::volScalarField & Zmax,
    Foam::volScalarField & Rmin, Foam::volScalarField & Rmax,
    Foam::volScalarField & Phimin, Foam::volScalarField & Phimax)
{
    const Foam::fvMesh & mesh = Zmin.mesh();
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        //---------------------------------------------------------
        const Foam::point & Position = mesh.faceCentres()[facei];
        const Foam::scalar Z = Position.z();
        const Foam::scalar R = Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y());
        const Foam::scalar phi = Angle(Position);
        //---------------------------------------------------------
        if(Zmin[owner[facei]] > Z) Zmin[owner[facei]] = Z;
        if(Zmax[owner[facei]] < Z) Zmax[owner[facei]] = Z;
        if(Rmin[owner[facei]] > R) Rmin[owner[facei]] = R;
        if(Rmax[owner[facei]] < R) Rmax[owner[facei]] = R;
        if(Phimin[owner[facei]] > phi) Phimin[owner[facei]] = phi;
        if(Phimax[owner[facei]] < phi) Phimax[owner[facei]] = phi;
        //---------------------------------------------------------
        if(Zmin[neighbour[facei]] > Z) Zmin[neighbour[facei]] = Z;
        if(Zmax[neighbour[facei]] < Z) Zmax[neighbour[facei]] = Z;
        if(Rmin[neighbour[facei]] > R) Rmin[neighbour[facei]] = R;
        if(Rmax[neighbour[facei]] < R) Rmax[neighbour[facei]] = R;
        if(Phimin[neighbour[facei]] > phi) Phimin[neighbour[facei]] = phi;
        if(Phimax[neighbour[facei]] < phi) Phimax[neighbour[facei]] = phi;
        //---------------------------------------------------------
    }
    forAll(mesh.boundaryMesh(), patchi)
    {
        const Foam::labelUList & pFaceCells = mesh.boundary()[patchi].faceCells();
        const Foam::vectorField & faceCentres = mesh.Cf().boundaryField()[patchi];
        if(faceCentres.empty()) continue;
        forAll(mesh.boundaryMesh()[patchi], facei)
        {
            //---------------------------------------------------------
            const Foam::vector & Position = faceCentres[facei];
            const scalar Z = Position.z();
            const Foam::scalar R = Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y());
            const scalar phi = Angle(Position);
            //---------------------------------------------------------
            if(Zmin[pFaceCells[facei]] > Z) Zmin[pFaceCells[facei]] = Z;
            if(Zmax[pFaceCells[facei]] < Z) Zmax[pFaceCells[facei]] = Z;
            if(Rmin[pFaceCells[facei]] > R) Rmin[pFaceCells[facei]] = R;
            if(Rmax[pFaceCells[facei]] < R) Rmax[pFaceCells[facei]] = R;
            if(Phimin[pFaceCells[facei]] > phi) Phimin[pFaceCells[facei]] = phi;
            if(Phimax[pFaceCells[facei]] < phi) Phimax[pFaceCells[facei]] = phi;
        }
    }
}

void UpdateContactField(const Foam::fvMesh & mesh,
    const Foam::volScalarField & Zmin, const Foam::volScalarField & Zmax,
    const Foam::volScalarField & Rmin, const Foam::volScalarField & Rmax,
    const Foam::volScalarField & Phimin, const Foam::volScalarField & Phimax,
    Foam::volScalarField & FieldZmin1, Foam::volScalarField & FieldZmax1, Foam::volScalarField & FieldZmin2, Foam::volScalarField & FieldZmax2,
    Foam::volScalarField & FieldRmin1, Foam::volScalarField & FieldRmax1, Foam::volScalarField & FieldRmin2, Foam::volScalarField & FieldRmax2,
    Foam::volScalarField & FieldPhimin, Foam::volScalarField & FieldPhimax,
    const Foam::surfaceScalarField & Field1f, const Foam::surfaceScalarField & Field2f, const Foam::surfaceScalarField & FieldPhif)
{
    const Foam::scalar Epsilon = Foam::SMALL*10.0;
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        //---------------------------------------------------------
        const Foam::point & Position = mesh.faceCentres()[facei];
        const Foam::scalar Z = Position.z();
        const Foam::scalar R = Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y());
        const Foam::scalar phi = Angle(Position);
        //---------------------------------------------------------
        if(Foam::mag(Zmin[owner[facei]] - Z) < Epsilon) FieldZmin1[owner[facei]] = Field1f[facei];
        if(Foam::mag(Zmax[owner[facei]] - Z) < Epsilon) FieldZmax1[owner[facei]] = Field1f[facei];
        if(Foam::mag(Zmin[owner[facei]] - Z) < Epsilon) FieldZmin2[owner[facei]] = Field2f[facei];
        if(Foam::mag(Zmax[owner[facei]] - Z) < Epsilon) FieldZmax2[owner[facei]] = Field2f[facei];
        if(Foam::mag(Rmin[owner[facei]] - R) < Epsilon) FieldRmin1[owner[facei]] = Field1f[facei];
        if(Foam::mag(Rmax[owner[facei]] - R) < Epsilon) FieldRmax1[owner[facei]] = Field1f[facei];
        if(Foam::mag(Rmin[owner[facei]] - R) < Epsilon) FieldRmin2[owner[facei]] = Field2f[facei];
        if(Foam::mag(Rmax[owner[facei]] - R) < Epsilon) FieldRmax2[owner[facei]] = Field2f[facei];
        if(Foam::mag(Phimin[owner[facei]] - phi) < Epsilon) FieldPhimin[owner[facei]] = FieldPhif[facei];
        if(Foam::mag(Phimax[owner[facei]] - phi) < Epsilon) FieldPhimax[owner[facei]] = FieldPhif[facei];
        //---------------------------------------------------------
        if(Foam::mag(Zmin[neighbour[facei]] - Z) < Epsilon) FieldZmin1[neighbour[facei]] = Field1f[facei];
        if(Foam::mag(Zmax[neighbour[facei]] - Z) < Epsilon) FieldZmax1[neighbour[facei]] = Field1f[facei];
        if(Foam::mag(Zmin[neighbour[facei]] - Z) < Epsilon) FieldZmin2[neighbour[facei]] = Field2f[facei];
        if(Foam::mag(Zmax[neighbour[facei]] - Z) < Epsilon) FieldZmax2[neighbour[facei]] = Field2f[facei];
        if(Foam::mag(Rmin[neighbour[facei]] - R) < Epsilon) FieldRmin1[neighbour[facei]] = Field1f[facei];
        if(Foam::mag(Rmax[neighbour[facei]] - R) < Epsilon) FieldRmax1[neighbour[facei]] = Field1f[facei];
        if(Foam::mag(Rmin[neighbour[facei]] - R) < Epsilon) FieldRmin2[neighbour[facei]] = Field2f[facei];
        if(Foam::mag(Rmax[neighbour[facei]] - R) < Epsilon) FieldRmax2[neighbour[facei]] = Field2f[facei];
        if(Foam::mag(Phimin[neighbour[facei]] - phi) < Epsilon) FieldPhimin[neighbour[facei]] = FieldPhif[facei];
        if(Foam::mag(Phimax[neighbour[facei]] - phi) < Epsilon) FieldPhimax[neighbour[facei]] = FieldPhif[facei];
        //---------------------------------------------------------
    }
    forAll(mesh.boundaryMesh(), patchi)
    {
        const Foam::labelUList & pFaceCells = mesh.boundary()[patchi].faceCells();
        const Foam::vectorField & faceCentres = mesh.Cf().boundaryField()[patchi];
        if(faceCentres.empty()) continue;
        const Foam::fvsPatchField<Foam::scalar> & BoundField1f = Field1f.boundaryField()[patchi];
        const Foam::fvsPatchField<Foam::scalar> & BoundField2f = Field2f.boundaryField()[patchi];
        const Foam::fvsPatchField<Foam::scalar> & BoundFieldPhif = FieldPhif.boundaryField()[patchi];
        forAll(mesh.boundaryMesh()[patchi], facei)
        {
            //---------------------------------------------------------
            const Foam::vector & Position = faceCentres[facei];
            const Foam::scalar Z = Position.z();
            const Foam::scalar R = Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y());
            const Foam::scalar phi = Angle(Position);
            //---------------------------------------------------------
            if(Foam::mag(Zmin[pFaceCells[facei]] - Z) < Epsilon) FieldZmin1[pFaceCells[facei]] = BoundField1f[facei];
            if(Foam::mag(Zmax[pFaceCells[facei]] - Z) < Epsilon) FieldZmax1[pFaceCells[facei]] = BoundField1f[facei];
            if(Foam::mag(Zmin[pFaceCells[facei]] - Z) < Epsilon) FieldZmin2[pFaceCells[facei]] = BoundField2f[facei];
            if(Foam::mag(Zmax[pFaceCells[facei]] - Z) < Epsilon) FieldZmax2[pFaceCells[facei]] = BoundField2f[facei];
            if(Foam::mag(Rmin[pFaceCells[facei]] - R) < Epsilon) FieldRmin1[pFaceCells[facei]] = BoundField1f[facei];
            if(Foam::mag(Rmax[pFaceCells[facei]] - R) < Epsilon) FieldRmax1[pFaceCells[facei]] = BoundField1f[facei];
            if(Foam::mag(Rmin[pFaceCells[facei]] - R) < Epsilon) FieldRmin2[pFaceCells[facei]] = BoundField2f[facei];
            if(Foam::mag(Rmax[pFaceCells[facei]] - R) < Epsilon) FieldRmax2[pFaceCells[facei]] = BoundField2f[facei];
            if(Foam::mag(Phimin[pFaceCells[facei]] - phi) < Epsilon) FieldPhimin[pFaceCells[facei]] = BoundFieldPhif[facei];
            if(Foam::mag(Phimax[pFaceCells[facei]] - phi) < Epsilon) FieldPhimax[pFaceCells[facei]] = BoundFieldPhif[facei];
        }
    }
}

inline bool IsContactZ(const Foam::scalar InductionZmin, const Foam::scalar InductionZmax, const Foam::scalar FieldZmin, const Foam::scalar FieldZmax, const Foam::scalar FieldPhimin, const Foam::scalar FieldPhimax, const Foam::vector U)
{
//    bool result = (FieldZmin < Foam::SMALL) && (FieldZmax > Foam::SMALL) && (Uphi < 120.0);    // Критерий через скорость
//    bool result = (FieldZmin < Foam::SMALL) && (FieldZmax > Foam::SMALL) && (FieldPhimax/FieldPhimin < 1.2); // Критерий через давление
    bool result = (InductionZmin < Foam::SMALL) && (InductionZmax > Foam::SMALL) && (FieldPhimax/FieldPhimin < 1.2); // Критерий через давление
    return result;
}


inline bool IsContactR(const Foam::scalar InductionRmin, const Foam::scalar InductionRmax, const Foam::scalar FieldRmin, const Foam::scalar FieldRmax, const Foam::scalar FieldPhimin, const Foam::scalar FieldPhimax, const Foam::vector U)
{
    bool result = (InductionRmin < Foam::SMALL) && (InductionRmax > Foam::SMALL) && (FieldPhimax/FieldPhimin < 1.2); // Критерий через давление
    return result;
}

void Contact(Foam::rhoRDEThermo & thermo, const Foam::volScalarField & Induction, Foam::volScalarField & MolWeight, const Foam::volVectorField & U, Foam::volScalarField & Upotential,
    const volScalarField & Zmin, const Foam::volScalarField & Zmax,
    const volScalarField & Rmin, const Foam::volScalarField & Rmax,
    const Foam::volScalarField & InductionZmin, const Foam::volScalarField & InductionZmax,
    const Foam::volScalarField & InductionRmin, const Foam::volScalarField & InductionRmax,
    const Foam::volScalarField & FieldZmin, const Foam::volScalarField & FieldZmax,
    const Foam::volScalarField & FieldRmin, const Foam::volScalarField & FieldRmax,
    const Foam::volScalarField & FieldPhimin, const Foam::volScalarField & FieldPhimax,
    Foam::scalarField & ContactCoordinateFieldZ, const Foam::vectorField & ContactCoordinateMeshZ,
    Foam::scalarField & ContactCoordinateFieldR, const Foam::vectorField & ContactCoordinateMeshR,
    const Foam::scalar dt, const Foam::scalar CameraLength, const Foam::scalar CameraDiameter)
{
    const Foam::fvMesh & mesh = MolWeight.mesh();
    // ПО ОСИ Z
    if(!ContactCoordinateMeshZ.empty())
        forAll(MolWeight, i) //if(Zmax[i] < CameraLength + Foam::SMALL) // если ячейка находится внутри камеры сгорания
            if(IsContactZ(InductionZmin[i], InductionZmax[i], FieldZmin[i], FieldZmax[i], FieldPhimin[i], FieldPhimax[i], U[i]))
    {
        // Получение индекса поля плоского "блинчика"(ContactCoordinateFieldZ), в котором хранятся значения контактных координат
        // index-i сетки камеры(mesh) => index-j сетки "блинчика"(ContactCoordinateMeshZ)
        const Foam::label j = GetLabelBounaryZ(mesh, ContactCoordinateMeshZ, i);
        // Получение ссылки на значение контактной координаты для i ячейки mesh и j ячейки ContactCoordinateFieldZ
        Foam::scalar & ContactCoordinate = ContactCoordinateFieldZ[j];
        //--------------------------------------------------------------------
        if(ContactCoordinate < Zmin[i]) ContactCoordinate = Zmin[i];
        if(ContactCoordinate > Zmax[i]) ContactCoordinate = Zmax[i];
        ContactCoordinate += dt*U[i].z();
        // alpha - доля остаточной горячей смеси, уменьшается до 0
        Foam::scalar alpha = (ContactCoordinate - Zmin[i]) / (Zmax[i] - Zmin[i]);
        if(alpha < 0.0) alpha = 0.0;
        if(alpha > 1.0) alpha = 1.0;
        // Интерполяция
        Foam::scalar Hthermodynamical = thermo.Hthermodynamical(i)*alpha + FieldZmax[i]*(1.0-alpha);
        Foam::ChemistryElement old = thermo.GetChemistry(i);
        Foam::scalar Hchemical = Upotential[i] - Hthermodynamical;
        // Кортеж из химических параметров
        ChemistryElement Element = thermo.GetChemistryFromHchemical(old, Hchemical);
        thermo.SetChemistry(i, Element);
//none        MolWeight[i] = MolWeight[i]*alpha + FieldZmax[i]*(1.0-alpha);
    }
//return;
    // ПО ОСИ R
    if(!ContactCoordinateFieldR.empty())
        forAll(MolWeight, i) //if(Zmax[i] < CameraLength + Foam::SMALL) // если ячейка находится внутри камеры сгорания
            if(IsContactR(InductionRmin[i], InductionRmax[i], FieldRmin[i], FieldRmax[i], FieldPhimin[i], FieldPhimax[i], U[i]))
    {
        // Получение индекса обертки(ContactCoordinateFieldR), в котором хранятся значения контактных координат
        // index-i сетки камеры(mesh) => index-j сетки "обертки"(ContactCoordinateMeshR)
        const Foam::label j = GetLabelBounaryR(mesh, ContactCoordinateMeshR, i);
        // Получение ссылки на значение контактной координаты для i ячейки mesh и j ячейки ContactCoordinateFieldR
        Foam::scalar & ContactCoordinate = ContactCoordinateFieldR[j];
        //--------------------------------------------------------------------
        if(ContactCoordinate < Rmin[i]) ContactCoordinate = Rmin[i];
        if(ContactCoordinate > Rmax[i]) ContactCoordinate = Rmax[i];
        const Foam::vector Centre = mesh.C()[i]; // Координата центра ячейки
        const Foam::scalar Ur = (U[i].x()*Centre.x() + U[i].y()*Centre.y()) / Foam::sqrt(Centre.x()*Centre.x() + Centre.y()*Centre.y());
        ContactCoordinate += dt*Ur;
        // alpha - доля остаточной горячей смеси, уменьшается до 0
        Foam::scalar alpha = (ContactCoordinate - Rmin[i]) / (Rmax[i] - Rmin[i]);
        if(alpha < 0.0) alpha = 0.0;
        if(alpha > 1.0) alpha = 1.0;
        // Интерполяция
        Foam::scalar Hthermodynamical = thermo.Hthermodynamical(i)*alpha + FieldRmax[i]*(1.0-alpha);
        Foam::ChemistryElement old = thermo.GetChemistry(i);
        Foam::scalar Hchemical = Upotential[i] - Hthermodynamical;
        // Кортеж из химических параметров
        ChemistryElement Element = thermo.GetChemistryFromHchemical(old, Hchemical);
        thermo.SetChemistry(i, Element);
//none        MolWeight[i] = MolWeight[i]*alpha + FieldRmax[i]*(1.0-alpha);
    }
}

inline bool CheckPointsRadius(const Foam::scalar CameraRadius, const Foam::pointField & Points)
{
    const Foam::scalar Epsilon = 1.0e-7;
    forAll(Points, i)
    {
        const Foam::scalar PointRadius = Foam::sqrt(Points[i].x()*Points[i].x() + Points[i].y()*Points[i].y());
        if(Foam::mag(CameraRadius - PointRadius) > Epsilon) return false;
    }
    return true;
}
#include <vector>
#include <utility>
// first - Z
// second - Mesh
typedef std::pair<Foam::scalar, Foam::vectorField> CellCandidate;
//-----------------------------------------------------------------------------------------------
// Функции для заполнения сетки(ContactCoordinateMesh) и данные(ContactCoordinateField) для контактной границы
// bool EmptyContactMesh = true, если в этой сетке нет области камеры сгорания
//-----------------------------------------------------------------------------------------------


void GetContactEdgeSize(const Foam::fvMesh & mesh, const Foam::scalar CameraLength, const Foam::scalar CameraDiameter, bool & EmptyContactMesh,
Foam::vectorField & ContactCoordinateMeshZ, Foam::scalarField & ContactCoordinateFieldZ,
Foam::vectorField & ContactCoordinateMeshR, Foam::scalarField & ContactCoordinateFieldR)
{
    const Foam::scalar CameraRadius = CameraDiameter / 2.0;
    // Z: Всю сетку режет на блины(или полоски) при Z=const, а затем выбирают блин с максимальным размером
    std::vector<CellCandidate> DATA;
    // резка на блины
    forAll(mesh.Cf(), i)
    {
        Foam::vector V = mesh.Cf()[i];
        if(V.z() > CameraLength + Foam::SMALL) continue; // Нет смысла обследовать грани, которые находятся выше камеры сгорания(диффузор-коллекторные)
        bool added = false;
        for(CellCandidate & elem : DATA) if(Foam::mag(V.z() - elem.first) < Foam::SMALL)
        {
            // Увеличиваем размер блина при определенным значения Z
            elem.second.append(V);
            added = true;
            break;
        }
        // Создаем новый блин
        if(!added)
        {
            CellCandidate cell;
            cell.first = V.z();
            cell.second.append(V);
            DATA.push_back(cell);
        }
    }
    // Если в этой сетке нет области камеры сгорания - чистим и выходим
    ContactCoordinateFieldZ.clear();
    if(DATA.empty())
    {
        EmptyContactMesh = true;
        return;
    }
    //----------------------------------------------
    // Выбираем блин с максимальным размером
    size_t CurrentZ = 0;
    for(size_t i = 0; i < DATA.size(); i++)
    if(DATA[CurrentZ].second.size() < DATA[i].second.size())
    {
        CurrentZ = i;
    }
    // Сетка с координатами
    ContactCoordinateMeshZ = DATA[CurrentZ].second;
    // Сетка с данными
    ContactCoordinateFieldZ.resize(ContactCoordinateMeshZ.size());
    forAll(ContactCoordinateFieldZ, i) ContactCoordinateFieldZ[i] = CameraLength;
    DATA.clear();

    
    //===================================================================================================
    // R: Всю сетку проверяем, если все точки(4шт) i-ой грани имеют значение R=CameraRadius, то кладем в ContactCoordinateMeshR
    ContactCoordinateMeshR.clear();
    forAll(mesh.faces(), i)
    {
        const Foam::point & FacePoint = mesh.faces()[i].centre(mesh.points());
        // Нет смысла обследовать грани, которые находятся выше камеры сгорания(диффузор-коллекторные)
        if(FacePoint.z() > CameraLength) continue;
        //-----------------------------------------------------
        // Получаемы поле точек iой грани. Всего должно быть 4 точки
        Foam::pointField PointFields = mesh.faces()[i].points(mesh.points());
        // Проверяем, сидят ли все точки на внешнем периметре
        // Если да, то грань - внешняя, её и запысываем в сетку
        if(CheckPointsRadius(CameraRadius, PointFields))
        {
            ContactCoordinateMeshR.append(FacePoint); // Добавим грань, у которой R=CameraRadius
        }
    }
    // Сетка с данными
    ContactCoordinateFieldR.resize(ContactCoordinateMeshR.size());
    forAll(ContactCoordinateFieldR, i) ContactCoordinateFieldR[i] = CameraRadius;
}
