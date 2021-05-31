/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(natural, [esyscmd(perl -e 'printf("%d", $1)')])
define(real, [esyscmd(perl -e 'printf("%0.2f", $1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


meshGenApp blockMesh;
convertToMeters 0.001;


// Настройки
define(Rtr, 1.0)// Коэфициент перед периодом
define(delta, 22)// Толщина зазора
define(diameter, 504)// Внешний диаметр
define(Length, 500)// Длина цилиндра
//----------------------------------------------------------
define(Rlarge, calc(Rtr * diameter / 2))// Максимальный радиус
define(Rsmall, calc(Rlarge - delta))// Минимальный радиус
define(NPerimeter, 400)// Количество ячееу по периметру
define(NRadius, 5)// Количество ячееу вдоль радиуса
define(NLength, 120)// Количество ячееу вдоль оси


// Конец настроек

define(PI, 3.1415926535897)
define(PositionXsmall, calc(Rsmall*cos(PI/4)))
define(PositionYsmall, calc(Rsmall*sin(PI/4)))
define(PositionXlarge, calc(Rlarge*cos(PI/4)))
define(PositionYlarge, calc(Rlarge*sin(PI/4)))
define(NSection, calc(NPerimeter/4))
define(Perimeter, calc(2.0*PI*Rlarge))
define(NN, natural(Perimeter))

/*
blockMesh параметры:
Минимальный радиус - Rsmall мм
Максимальный радиус - Rlarge мм
Длина цилиндра - Length мм
Внешний периметр цилиндра - Perimeter мм
Количество ячееу по периметру - NPerimeter
Количество ячееу вдоль радиуса - NRadius
Количество ячееу вдоль оси - NLength
*/
vertices
(
    ( Rsmall  0.0 0.0) // 0
    ( 0.0  Rsmall 0.0) // 1
    (-Rsmall  0.0 0.0) // 2
    ( 0.0 -Rsmall 0.0) // 3

    ( Rlarge  0.0 0.0) // 4
    ( 0.0  Rlarge 0.0) // 5
    (-Rlarge  0.0 0.0) // 6
    ( 0.0 -Rlarge 0.0) // 7

    ( Rsmall  0.0 Length) // 8
    ( 0.0  Rsmall Length) // 9
    (-Rsmall  0.0 Length) // 10
    ( 0.0 -Rsmall Length) // 11

    ( Rlarge  0.0 Length) // 12
    ( 0.0  Rlarge Length) // 13
    (-Rlarge  0.0 Length) // 14
    ( 0.0 -Rlarge Length) // 15
);
blocks
(
    // I четверть - верх-право
    hex (0 4 5 1 8 12 13 9) blockUPRIGHT (natural(NRadius) natural(NSection) natural(NLength)) simpleGrading (1 1 1)
    // II четверть - верх-лево
    hex (1 5 6 2 9 13 14 10) blockUPLEFT (natural(NRadius) natural(NSection) natural(NLength)) simpleGrading (1 1 1)
    // III четверть - низ-лево
    hex (2 6 7 3 10 14 15 11) blockDOWNLEFT (natural(NRadius) natural(NSection) natural(NLength)) simpleGrading (1 1 1)
    // IV четверть - низ-право
    hex (3 7 4 0 11 15 12 8) blockDOWNRIGHT (natural(NRadius) natural(NSection) natural(NLength)) simpleGrading (1 1 1)
);

//create the quarter circles
edges
(
    // Дуги для Z = 0
    // I четверть - верх-право
    arc 0 1 ( real(PositionXsmall)  real(PositionYsmall) 0.0)
    arc 4 5 ( real(PositionXlarge)  real(PositionYlarge) 0.0)
    // II четверть - верх-лево
    arc 1 2 ( real(-PositionXsmall)  real(PositionYsmall) 0.0)
    arc 5 6 ( real(-PositionXlarge)  real(PositionYlarge) 0.0)
    // III четверть - низ-лево
    arc 2 3 ( real(-PositionXsmall) real(-PositionYsmall) 0.0)
    arc 6 7 ( real(-PositionXlarge) real(-PositionYlarge) 0.0)
    // IV четверть - низ-право
    arc 3 0 ( real(PositionXsmall) real(-PositionYsmall) 0.0)
    arc 7 4 ( real(PositionXlarge) real(-PositionYlarge) 0.0)

    // Дуги для Z = Length
    // I четверть - верх-право
    arc 8 9  ( real(PositionXsmall)  real(PositionYsmall) Length)
    arc 12 13( real(PositionXlarge)  real(PositionYlarge) Length)
    // II четверть - верх-лево
    arc 9 10 ( real(-PositionXsmall)  real(PositionYsmall) Length)
    arc 13 14( real(-PositionXlarge)  real(PositionYlarge) Length)
    // III четверть - низ-лево
    arc 10 11( real(-PositionXsmall) real(-PositionYsmall) Length)
    arc 14 15( real(-PositionXlarge) real(-PositionYlarge) Length)
    // IV четверть - низ-право
    arc 11 8 ( real(PositionXsmall) real(-PositionYsmall) Length)
    arc 15 12( real(PositionXlarge) real(-PositionYlarge) Length)
);

boundary
(
    inlet // Когда Z = Length
    {
        type patch;
        faces
        (
            (8 12 13 9)
            (9 13 14 10)
            (10 14 15 11)
            (11 15 12 8)
        );
    }
    outlet // Когда Z = 0
    {
        type patch;
        faces
        (
            (0 4 5 1)
            (1 5 6 2)
            (2 6 7 3)
            (3 7 4 0)
        );
    }
    walls
    {
        type wall;
        faces
        (
            // I четверть - верх-право
            (0 1 9 8)
            (4 5 13 12)
            // II четверть - верх-лево
            (1 2 10 9)
            (5 6 14 13)
            // III четверть - низ-лево
            (2 3 11 10)
            (6 7 15 14)
            // IV четверть - низ-право
            (3 0 8 11)
            (7 4 12 15)
        );
    }
);
mergePatchPairs
(
);