// Список скалярных переменных
// 0 - время сущствования жесткой стенки
// 1 - Backpressure
// 2 - TimeKoefficientDetonationPerimeter
// 3 - CameraLength
// 4 - CameraDiameter
// торец
// 5 - SstarS
// 6 - Pstar
// 7 - Tstar
// боковая стенка
// 8 - SstarS
// 9 - Pstar
// 10 - Tstar
// 11 - min deltaT
Foam::scalarIOList scalarParameters
(
    IOobject
    (
        "scalarParameters",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    12
);


Foam::dictionary dictionaryBoundary("./constant/Properties");
std::filebuf fb;
bool b = fb.open(dictionaryBoundary.name(), std::ios::in);
std::istream stream(&fb);
Foam::ISstream Sstream(stream, dictionaryBoundary.name());
dictionaryBoundary.read(Sstream);
// Чтение
scalarParameters[0] = 0.0;
scalarParameters[1] = readScalar(dictionaryBoundary.subDict("geometry").lookup("Backpressure"));
scalarParameters[2] = readScalar(dictionaryBoundary.subDict("geometry").lookup("TimeKoefficientDetonationPerimeter"));
scalarParameters[3] = readScalar(dictionaryBoundary.subDict("geometry").lookup("CameraLength"));
scalarParameters[4] = readScalar(dictionaryBoundary.subDict("geometry").lookup("CameraDiameter"));
// торец
scalarParameters[5] = readScalar(dictionaryBoundary.subDict("butt").lookup("SstarS"));
scalarParameters[6] = readScalar(dictionaryBoundary.subDict("butt").lookup("Pstar"));
scalarParameters[7] = readScalar(dictionaryBoundary.subDict("butt").lookup("Tstar"));
// боковая стенка
scalarParameters[8] = readScalar(dictionaryBoundary.subDict("sidewall").lookup("SstarS"));
scalarParameters[9] = readScalar(dictionaryBoundary.subDict("sidewall").lookup("Pstar"));
scalarParameters[10] = readScalar(dictionaryBoundary.subDict("sidewall").lookup("Tstar"));
// Прочее
Tmin = readScalar(dictionaryBoundary.subDict("geometry").lookup("Tmin"));
scalarParameters[11] = readScalar(dictionaryBoundary.subDict("geometry").lookup("deltaTmin"));
// Не вещественные параметры
//----------------------------------
Multiplier = readLabel(dictionaryBoundary.subDict("geometry").lookup("Multiplier"));        // Multiplier - умножитель волн
bool UseChemistry = readBool(dictionaryBoundary.subDict("geometry").lookup("UseChemistry"));   // UseChemistry - подключать ли химию
bool UseContact = readBool(dictionaryBoundary.subDict("geometry").lookup("UseContact"));       // Использовать ли контактную границу
ReadBreakDown(dictionaryBoundary.subDict("geometry"));
// Определение THE - Temperature-Hentalpy-Energy
Foam::word StringTHE = dictionaryBoundary.subDict("geometry").lookup("Tuplayer");
if(StringTHE == "THE") UpLayer = EUpLayer::THE;
else if(StringTHE == "getT") UpLayer = EUpLayer::getT;
else UpLayer = EUpLayer::THE;
// Закрытие
fb.close();

// Время, за которое детонационная волна пройдет один периметр
Foam::scalar & timeOnePerimeter = scalarParameters[0];
// Противодавление
dimensionedScalar Backpressure("Backpressure", dimensionSet(1, -1, -2, 0, 0, 0, 0), scalarParameters[1]);
// Коэфициент перед временем, которое необходимо для прохождения ДВ со скоростью Ч-Ж за переметр
const Foam::scalar & TimeKoefficientDetonationPerimeter = scalarParameters[2];
// Длина камеры сгорания
CameraLength = scalarParameters[3];
// Диаметр камеры сгорания
CameraDiameter = scalarParameters[4];
// Минимальный шаг по времени
const Foam::scalar & deltaTmin = scalarParameters[11];


