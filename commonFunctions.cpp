#include "commonFunctions.h"
//------------------------------------------------------------------------------
inline MyDataType dist(MyDataType x1, MyDataType x2)
{
    return fabsl(x1-x2);
}

void calculatePolinomByKoef(TSignal & inX, TSignal & koef, TSignal & outF )
{
    auto NumberOfPoints=inX.size();
    auto NumberOfPolinomKoef=koef.size()-1;

    outF.clear();

    for( auto i=0u;i<NumberOfPoints;++i)
    {
        outF.push_back(0);
        for (auto koef_index = 0u; koef_index <= NumberOfPolinomKoef; ++koef_index)
        {
            long double powedB=0;
            if(NumberOfPolinomKoef-koef_index==0)
            powedB=1;
            else
            powedB=powl(inX[i],NumberOfPolinomKoef-koef_index);

            outF[i]+=koef[koef_index]*powedB;
        }
    }

}

void SaveToFile(std::vector<std::string> *tsl, std::string filename)
{
    std::ofstream fout(filename, std::ofstream::out);

    for(auto s=tsl->begin();s!=tsl->end();++s)
    {
        fout << *s<< "\n";
    }
    fout.close();
}

bool LoadFromFile(std::vector<std::string> *tsl, std::string filename)
{
    std::string t;
    std::ifstream fin(filename);
    while(fin.good())
    {
        fin >> t;
        tsl->push_back(t);
    }
    return true;
}


bool testCommonFunctions()
{
    long double e=10e-6L;
    long double x=1.123456L;
    long double y=1.123454L;
    unsigned int NumberOfDecimalPlaces=5;
    if(fabsl(roundM(x,NumberOfDecimalPlaces)-1.12346L)>e)
    {
        return false;
    }
    if (fabsl(roundM(y,NumberOfDecimalPlaces)-1.12345L)>e)
    {
        return false;
    }

    TSignal s;
    s.push_back(1.123456L);
    s.push_back(1.123454L);
    roundM(s,NumberOfDecimalPlaces);
    if (fabsl(s[0]-1.12346L)>e || fabsl(s[1]-1.12345L)>e)
    {
        return false;
    }
    return true;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
TSignal calculateSecondDerivative(TSignal & y, MyDataType h)
{
  TSignal d2Y(y.size()-2);

      for(size_t i =0;i<y.size()-2;i++)
      {
        //d2Y[i]=1.0/(fabs(dY[i+1]-dY[i]))*(2.0*dY[i+1]-dY[i+2]/2.0-3.0/2.0*dY[i]);
        // формула f(x-h)-2f(x)+f(x+h)/h^2
        d2Y[i]=(y[i]-2*y[i+1]+y[i+2])/h/h;
        //d2Y[i]=(dY[i+2]-dY[i])/2.0/h;
      }
  return d2Y;
}

//-------------------------------------------------------------------
TSignal calculateFirstDerivative(TSignal & y, MyDataType h)
{
  TSignal dY(y.size()-2);
  // Посчитаем производную методом конечных разностей
  // формула df/dx=(f(x+h)-f(x-h))/2/h;
  for(size_t i =0;i<y.size()-2;i++)
      {
        //dY[i]=1.0/(fabs(y[i+1]-y[i]))*(2.0*y[i+1]-y[i+2]/2.0-3.0/2.0*y[i]);
        dY[i]=(y[i+2]-y[i])/2.0L/h;
      }
  return dY;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

long double calculateMeanValue(const TSignal & x)
{
    long double res=0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        res+=x[i];
    }
    return res/x.size();
}

long double calculateDispersion(const TSignal & x)
{
    long double res=0;
    long double mean=calculateMeanValue(x);
    for (size_t i = 0; i < x.size(); ++i)
    {
        res+=powl(x[0]-mean,2);
    }
    return res;
}


long double standardDeviation(const TSignal & x)
{
    return sqrtl(calculateDispersion(x)/x.size());
}

//------------------------------------------------------------------------------



//---------------------------------------------------------------------------
// округление результатов при сохранении, необходимо потому что если этого не сделать
// метод будет считать что мы с этой точностью измеряли - и мало что найдет
void RoundM(long double * x,int length, size_t NumberOfNumbersAfterPoint)
{
    int S=powl(10,NumberOfNumbersAfterPoint);
    for(int i=0;i<length;i++)
    {
        int n=static_cast<int>(x[i]*S)%10;
        if(n<5)
            x[i]=static_cast<long double>(floorl(x[i]*S))/S;
        else
            x[i]=static_cast<long double>(ceill(x[i]*S))/S;
    }
}
// Округление с заданной точностью.
void roundM(TSignal & in, unsigned int & NumberOfDecimalPlaces )
{

    long long S=pow(10,NumberOfDecimalPlaces);
    for(auto pos=in.begin();pos!=in.end();++pos)
    {
        int n=static_cast<int>(*pos*S*10)%10;
        if(n<5.0)
            *pos=static_cast<long double>(floorl(*pos*S))/S;
        else
            *pos=static_cast<long double>(ceill(*pos*S))/S;
    }
}
//-------------------------------------------------------------------------------

MyDataType roundM(MyDataType x, unsigned int NumberOfDecimalPlaces )
{
    long double S=pow(10,NumberOfDecimalPlaces);
        MyDataType n=static_cast<int>(x*S*10.0)%10;
        if(n<5.0)
            x=static_cast<long double>(floorl(x*S))/S;
        else
            x=static_cast<long double>(ceill(x*S))/S;
    return x;
}
//-------------------------------------------------------------------------------


//------------------------------------------------------------------------------
/*
Прореживание зависимостей - по идее функция универсальная.
Способная из массива одной длины получить массив меньшей длины.
При этом рассчитывается новый равномерный шаг по сетке.

Используется для:
1. Согласования измеренных данных после получения их с АЦП.
2. Согласования данных после фильтрации.
3. Согласования данных после экстраполяции.
4. Сохранении данных с заданным шагом и количеством точек.


Есть такая проблема - сигнал для отрицательного магнитного поля поступает сюда
в порядке от 0 до -2
т.е. вообще говоря в обратном, поэтому это не пашет.
*/
bool thiningSignal(TSignal & inB, TSignal & inDependence, TSignal & outB, TSignal & outDependence,
    MyDataType left, MyDataType right, size_t NewLength)
{
    if (right<left) // если при вызове перепутали границы.
    {
        MyDataType temp=right;
        right = left;
        left = temp;
    }

    size_t OldLength = inB.size();

    if (NewLength>OldLength) // если просят увеличить количество точек.
    {
        return false;
    }

    if (OldLength != inDependence.size())
    {
        return false; // если зависимости разной длины - мы не можем им помочь.
    }

    if (OldLength==NewLength)
    {
    // при совпадении размеров - нам делать совсем нечего.
        outDependence=inDependence;
        outB=inB;
        return true;
    }

    outB.clear(); // чистим выходные зависимости, на случай если это не сделали за нас.
    outDependence.clear();

    TSignal idealB; // тут будем хранить опорные точки.

    // шаг есть величина диапазона на количество интервалов (на единицу меньше количества точек)
    MyDataType shag=(right-left)/(static_cast<MyDataType>(NewLength)-1.0L);

    idealB.push_back(left); // начинаем с наименьшей границы
    for (unsigned int i=1; i < NewLength; ++i) 
    {
    // создаем сетку для магнитного поля.
        idealB.push_back(idealB[i-1]+shag);
    }
    // тут начинается поиск.
    /*
    TSignal indexes(NewLength);
    TSignal distance(NewLength);
    for (int i = 0; i < NewLength; ++i)
    {
        indexes[i]=0;
        distance[i]=dist(inB[i],idealB[i]);
    }

    for (int i = 1; i < NewLength; ++i)
    {
        if (dist(inB[i],idealB[i])<distance[i])
        {
            distance[i]=dist(inB[i],idealB[i]);
            indexes[i]=i;
        }
    }*/
    // возможно стоит прикрутить более оптимальный поиск.
    /*
    for (unsigned int i = 0; i < NewLength; ++i)
    {
    // одно из самых тупых мест в этой программе.
    // может прикрутить сюда бинарный поиск?
    // данные-то по идее отсортированы
    // не считая точек возле нулевого значения магнитной индукции.
        unsigned int index=0;
        long double r=4;
        for(unsigned int k=0;k<OldLength;++k)
        {
            if(dist(inB[k],idealB[i])<=r)
            {
                r=dist(inB[k],idealB[i]);
                index=k;
            }
        }
        if(index<OldLength)
        {
            outB.push_back(inB[index]);
            outDependence.push_back(inDependence[index]);
        }
    }
       */
    for (unsigned int i = 0; i < NewLength; ++i)
    {

    if(inB[0]>inB[inB.size()/2]) // попробуем отделить отрицательный сигнал.
    {    // нас волнует именно порядок элементов
    // т.е. в начале у нас в районе нуля
    // в конце в районе -2
    std::reverse(inB.begin(),inB.end());
    }
    std::vector<MyDataType>::iterator b=std::lower_bound(inB.begin(),inB.end(),idealB[i]);
    if(b==inB.end())
        --b;
    if(b!=inB.begin())
    {
        if (fabsl(fabsl(*b)-fabsl(idealB[i]))>fabsl(fabsl(*(b-1))-fabsl(idealB[i])))
            --b;
    }
    outB.push_back(*b);
    unsigned long d=static_cast<unsigned long>(std::distance(inB.begin(),b));
    outDependence.push_back(inDependence[d]);

    }

    return true;
}
//------------------------------------------------------------------------------
