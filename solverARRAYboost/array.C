#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "fvCFD.H"
#include "primitiveFields.H"

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
using namespace boost::filesystem;

// Структура ячейки: время - положение пика - значение пика
struct SCell
{
    Foam::scalar Time;
    Foam::scalar Xpeak;
    Foam::label ipeak;
    Foam::scalar peak;
    Foam::scalar D;
    std::string PathName;
    // CJ
    Foam::scalar XCJ;
    Foam::scalar c;
    Foam::scalar U;
    Foam::scalar T;
    Foam::scalar p;
    Foam::scalar rho;
    Foam::scalar molWeight;
};
// База данных ячеек
std::vector<SCell> DATA;

void MAX(const Foam::fvMesh & mesh, const Foam::scalarField & Field, SCell & Cell)
{
    Cell.peak = 0.0;
    forAll(Field, i) if(Field[i] > Cell.peak)
    {
        Cell.peak = Field[i];
        Cell.Xpeak = mesh.C()[i].x();
        Cell.ipeak = i;
    }
}
inline void ReadLines(std::istream & stream, size_t count)
{
    while(count > 0)
    {
        char ch = stream.get();
        if(ch == '\n') count--;
    }
}

void directory_use(const Foam::fvMesh & mesh, std::string DirectoryName, std::string FileName, std::string PathName)
{
    SCell Cell;
    Cell.Time = std::atof(DirectoryName.c_str());
    Cell.PathName = PathName;
    
    std::string file = PathName + "/" + FileName;
    // Читалка для поля
    std::filebuf FileBuffer;
    bool b = FileBuffer.open(file, std::ios::in);
    std::istream stream(&FileBuffer);
    Foam::ISstream Sstream(stream, file);
    ReadLines(Sstream.stdStream(), 20);
    // Само поле
    Foam::scalarField Field(Sstream);
    FileBuffer.close();
    MAX(mesh, Field, Cell);
    DATA.push_back(Cell);
}

void CJ(const Foam::fvMesh & mesh, SCell & Cell)
{
    std::filebuf buffer;
    // U
    buffer.open(Cell.PathName + "/U", std::ios::in);
    std::istream Ustream(&buffer);
    Foam::ISstream USstream(Ustream, Cell.PathName + "/U");
    ReadLines(USstream.stdStream(), 20);
    Foam::vectorField U(USstream); // Само поле
    buffer.close();
    // c
    buffer.open(Cell.PathName + "/c", std::ios::in);
    std::istream Cstream(&buffer);
    Foam::ISstream CSstream(Cstream, Cell.PathName + "/c");
    ReadLines(CSstream.stdStream(), 20);
    Foam::scalarField c(CSstream); // Само поле
    buffer.close();
    // T
    buffer.open(Cell.PathName + "/T", std::ios::in);
    std::istream Tstream(&buffer);
    Foam::ISstream TSstream(Tstream, Cell.PathName + "/T");
    ReadLines(TSstream.stdStream(), 20);
    Foam::scalarField T(TSstream); // Само поле
    buffer.close();
    // p
    buffer.open(Cell.PathName + "/p", std::ios::in);
    std::istream Pstream(&buffer);
    Foam::ISstream PSstream(Pstream, Cell.PathName + "/p");
    ReadLines(PSstream.stdStream(), 20);
    Foam::scalarField p(PSstream); // Само поле
    buffer.close();
    // rho
    buffer.open(Cell.PathName + "/rho", std::ios::in);
    std::istream Rhostream(&buffer);
    Foam::ISstream RhoSstream(Pstream, Cell.PathName + "/rho");
    ReadLines(RhoSstream.stdStream(), 20);
    Foam::scalarField rho(PSstream); // Само поле
    buffer.close();
    // molWeight
    buffer.open(Cell.PathName + "/MolWeight", std::ios::in);
    std::istream molWeightstream(&buffer);
    Foam::ISstream molWeightSstream(molWeightstream, Cell.PathName + "/MolWeight");
    ReadLines(molWeightSstream.stdStream(), 20);
    Foam::scalarField molWeight(molWeightSstream); // Само поле
    buffer.close();
    //-------------------------------------------------------------
    // Поиск Чепмена-Жуге
    Foam::label iCJ = 0;
    Foam::label imax = U.size()-1;
    /*Foam::scalar min = Foam::GREAT;
    forAll(U, i) if(Foam::mag(U[i].x() + c[i] - Cell.D) < min)
    {
        iCJ = i;
        min = Foam::mag(U[i].x() + c[i] - Cell.D);
    }*/
    //-------------------------------------------------------------
    bool DmoreUpC = false; // Курок сброса
    for(Foam::label i = imax; i >= 0; i--)
    {
        if(!DmoreUpC && (U[i].x() + c[i] > Cell.D))
        {
            // Взводим курок
            DmoreUpC = true;
            continue;
        }
        if(DmoreUpC && (U[i].x() + c[i] < Cell.D))
        {
            // Спускаем курок
            iCJ = i;
            break;
        }
    }
    Cell.XCJ = mesh.C()[iCJ].x();
    Cell.U = U[iCJ].x();
    Cell.T = T[iCJ];
    Cell.c = c[iCJ];
    Cell.p = p[iCJ];
    Cell.rho = rho[iCJ];
    Cell.molWeight = molWeight[iCJ];
}

bool predicate(const SCell & a, const SCell & b)
{
    return a.Time < b.Time;
}
int main(int argc, char *argv[])
{
    std::string FileName("p");
    if(argc > 1) FileName = argv[1];
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // Получаем теущий путь
    boost::filesystem::path thispath = boost::filesystem::current_path();
    // Цикл по всем файлам в текущей директории(причем поддиректории - это тоже файл, что нам и надо)
    for (directory_iterator file(thispath); file != directory_iterator(); ++file)
    {
//        if(is_directory(file)) // если iй файл - поддиректория
        {
    	    path p = *file;
            // Итак, знаем, что текущий файл - директория, но нужно проверить, является ли имя файла числом.
            std::string DirectoryName = p.string().substr(p.string().rfind("/") + 1, std::string::npos);
            // Если строку никак не преобразовать в этот тип, то возвращается 0.0
            double test = std::atof(DirectoryName.c_str());
            if(test != 0.0) directory_use(mesh, DirectoryName, FileName, p.string()); // Обработка папки
        }
    }
    // Сортируем с использованием лямда-выражения
//    std::sort(DATA.begin(), DATA.end(), [=](SCell & a, SCell & b)->bool{return a.Time < b.Time;});
    std::sort(DATA.begin(), DATA.end(), predicate);
    // Вычисляем D
    Foam::scalar told = 0.0;
    Foam::scalar Xold = 0.0;
    for(std::vector<SCell>::iterator iter = DATA.begin(); iter != DATA.end(); iter++)
    {
        if(iter == DATA.begin())
        {
            (*iter).D = ((*iter).Xpeak - 0.0) / ((*iter).Time - 0.0);
            told = (*iter).Time;
            Xold = (*iter).Xpeak;
        }else
        {
            (*iter).D = ((*iter).Xpeak - Xold) / ((*iter).Time - told);
            told = (*iter).Time;
            Xold = (*iter).Xpeak;
        }
    }
    // Вычисляем Чепмен-Жуге
    for(SCell & Cell : DATA) CJ(mesh, Cell);
    // Записываем результат
    std::ofstream out("data.dat");
    out << "#        time    -    Xpeak   -   peak(atm)-   D    -    XCJ   -    c       -      U    -    molWeight   -  T     -   p(pas)     -   p(atm)     -   rho" << std::endl;
    for(SCell & cell : DATA)
    {
        out.width(14); out << cell.Time;
        out.width(14); out << cell.Xpeak;
        out.width(14); out << cell.peak;
        out.width(8); out << cell.D;
        out.width(12); out << cell.XCJ;
        out.width(12); out << cell.c;
        out.width(14); out << cell.U;
        out.width(12); out << cell.molWeight;
        out.width(12); out << cell.T;
        out.width(14); out << cell.p;
        out.width(14); out << cell.p/101325.0;
        out.width(14); out << cell.rho;
        out << std::endl;
    }
    Info << "END WORDKS" << endl;
    return 0;
}
