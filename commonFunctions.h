#ifndef COMMONFUNCTIONSH
#define COMMONFUNCTIONSH


#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <set>
#include <iostream>
#include <fstream>

typedef long double MyDataType;
typedef std::vector<MyDataType> TSignal;
typedef TSignal::iterator TSignalIt;
typedef TSignal::const_iterator TSignalconstIt;
typedef std::vector<TSignal> TwoDimensionalContainer;
typedef std::vector<std::vector<MyDataType> > TMatrix;

typedef std::vector<long double> InDataSpectr;
typedef std::vector<std::vector<long double> > TStatistic;
typedef std::vector<long double> MyData_spektr;

typedef std::vector<std::string> TStringList;


enum SaveType {ALL_POINTS,POINTS_11,POINTS_21};

void SaveToFile(std::vector<std::string> *tsl, std::string filename);
bool LoadFromFile(std::vector<std::string> *tsl, std::string filename);

void RoundM(long double * x,int length, size_t NumberOfNumbersAfterPoint);
void roundM(TSignal & in, unsigned int & NumberOfDecimalPlaces );
MyDataType roundM(MyDataType x, unsigned int NumberOfDecimalPlaces );
void calculatePolinomByKoef(TSignal & inX, TSignal & koef, TSignal & outF );

TSignal calculateSecondDerivative(TSignal & y, MyDataType h);
TSignal calculateFirstDerivative(TSignal & y, MyDataType h);
long double calculateMeanValue(const TSignal & x);
long double calculateDispersion(const TSignal & x);
long double standardDeviation(const TSignal & x);

bool testCommonFunctions();

struct PeaksCriteria
{
  MyDataType peakHeigh;
  MyDataType peakWidth;
  MyDataType peakWidthOrd;
  MyDataType peakVelocityR;
  MyDataType peakVelocity2R;
  MyDataType symmetri;
  MyDataType peakVelocityL;
  MyDataType peakVelocity2L;
  MyDataType peakLeftBorderFirstDerivative;
  MyDataType peakFirstDerivative;
  MyDataType peakRightBorderFirstDerivative;
  MyDataType peakLeftMiddleFirstDerivative;
  MyDataType peakRightMiddleFirstDerivative;
  MyDataType peakLeftBorderSecondDerivative;
  MyDataType peakSecondDerivative;
  MyDataType peakRightBorderSecondDerivative;
  MyDataType peakLeftMiddleSecondDerivative;
  MyDataType peakRightMiddleSecondDerivative;
};

struct AdditionalData
{
  MyDataType leftPointElectronConductivity;
  MyDataType leftPointHoleConductivity;
  MyDataType rightPointElectronConductivity;
  MyDataType rightPointHoleConductivity;
  MyDataType leftPointElectronConductivityLog;
  MyDataType leftPointHoleConductivityLog;
  MyDataType rightPointElectronConductivityLog;
  MyDataType rightPointHoleConductivityLog;
  MyDataType leftPointMobility;
  MyDataType rightPointMobility;
};

/*
Давайте поиграем в продумывание архитектуры
Вот есть у нас завивимости:
Поле
и Какие-то сигналы
Поле в общем-то тоже сигнал

Как результат наш набор данных должен быть массивом сигналов?
С явными обозначениями, так что видимо это действительно будет
либо ассоциативный массив map либо набор пар

Получается мне нужно реализовать отдельно класс сигнала в принципе
а потом наследовать его и переопределять некоторые функции согласно тому, как должен обрабатываться конкретный сигнал...
я в шоке... профессиональная деформация началась.

*/
typedef std::set<std::pair<std::string, TSignal> >DataSet;

typedef std::vector< std::pair<std::string, std::string> > optionDescription;

enum DataKind {CURRENT_DATA=0, FILTERED_DATA=1, EXTRAPOLATED_DATA=2, AVERAGED_DATA=3, TENZOR_DATA=4, MOBILITY_DATA=5, FIT_DATA=6};

enum SignalType {MAGNETIC_FIELD=0, HALL_EFFECT=1, MAGNETORESISTANCE=2, PHOTOCONDUCTIVITY_FARADEY=3, PHOTOCONDICTIVITY_FOYGHT=4,MAGNETIC_FIELD_F=5, SXX=6,SXY=7};
enum ParamsType {DIRECT=0,REVERSE=1,COMBINE=2};

const long double THEALMOSTZERO = 0.000001;


inline MyDataType dist(MyDataType x1, MyDataType x2);

bool thiningSignal(TSignal & inB, TSignal & inDependence, TSignal & outB, TSignal & outDependence,
    MyDataType left, MyDataType right, size_t NewLength);


//---------------------------------------------------------------------------


template <class T>
size_t indexOfElemClosestToValue(const std::vector<T> & in,MyDataType value)
{
    MyDataType dist=fabs(in[0]-value);
    size_t index=0;
    for (size_t i = 1; i < in.size(); ++i)
    {
        if (dist>fabs(in[i]-value))
        {
            dist=fabs(in[i]-value);
            index=i;
        }
    }
    return index;
}
//-------------------------------------------------------------------------------

template <class T>
size_t indexOfElemClosestToValueFromLeft(const std::vector<T> & in,MyDataType value)
{
    MyDataType dist=fabs(in[0]-value);
    size_t index=0;
    for (int i = 1; i < in.size(); ++i)
    {
        if (dist>fabs(in[i]-value))
        {
            dist=fabs(in[i]-value);
            index=i;
        }
    }
    return index;
}
//-------------------------------------------------------------------------------

template <class T>
size_t indexOfElemClosestToValueFromRight(const std::vector<T> & in,MyDataType value)
{
    MyDataType dist=fabs(in[0]-value);
    size_t index=0;
    for (int i = in.size()-1; i >= 0; --i)
    {
        if (dist>fabs(in[i]-value))
        {
            dist=fabs(in[i]-value);
            index=i;
        }
    }
    return index;
}

//-------------------------------------------------------------------------------
template <class T>
inline T max_elem(const std::vector<T> & in)
{
    if(in.size()==0)
    {
        //"Сложно найти максимальный элемент, когда нет ни одного");
        return 0;
    }
    T max=*in.begin();
    for(auto pos=in.begin();pos!=in.end();++pos)
    {
        if(*pos>max)
        max=*pos;
    }
    return max;
}
//-------------------------------------------------------------------------------
template <class T>
inline T max_abs_elem(const std::vector<T> & in)
{
    if(in.size()==0)
    {
        //"Сложно найти максимальный элемент, когда нет ни одного");
        return 0;
    }
    T max=fabs(*in.begin());
    for(auto pos=in.begin();pos!=in.end();++pos)
    {
        if(fabs(*pos)>max)
        max=fabs(*pos);
    }
    return max;
}
//-------------------------------------------------------------------------------

template <class T>
inline T min_elem(const std::vector<T> & in)
{
    if(in.size()==0)
    {
        //"Сложно найти максимальный элемент, когда нет ни одного");
        return 0;
    }
    T min=*in.begin();
    for(auto pos=in.begin();pos!=in.end();++pos)
    {
        if(*pos<min)
        min=*pos;
    }
    return min;
}
//-------------------------------------------------------------------------------

template <class T>
inline T min_abs_elem(const std::vector<T> & in)
{
    if(in.size()==0)
    {
        //"Сложно найти максимальный элемент, когда нет ни одного");
        return 0;
    }
    T min=fabs(*in.begin());
    for(auto pos=in.begin();pos!=in.end();++pos)
    {
        if(fabs(*pos)<min)
        min=fabs(*pos);
    }
    return min;
}
//-------------------------------------------------------------------------------
template <class T>
int sign(T in)
{
    if (in>=0)
    {
        return 1;
    }
    else 
        return -1;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
Медианный фильтр.
Просто возвращает медиану массива.
*/
template <class T>
T medianFilter(std::vector <T> & in)
{
    std::sort(in.begin(),in.end());
    return in[in.size()/2];
}

template <class T>
T medianFilter(typename std::vector <T>::iterator &beg, typename std::vector <T>::iterator &end)
{
    std::sort(beg,end);
    return *(beg+std::distance(beg,end)/2);
}

template <class T>
void medianFilter(std::vector <T> & in,std::vector <T> & out,size_t size)
{
    if (in.size()<size)
    {
        out.push_back(medianFilter(in));
    }
    for(unsigned int i=size;i<in.size() ;i+=size)
    {
        std::sort(&in[i-size],&in[i]);
        out.push_back(in[i-size/2]);
    }
}

template <class T>
std::vector<T> calculateSecondDerivative(std::vector<T> & y, T h)
{
  std::vector<T> d2Y(y.size()-2);

      for(int i =0;i<y.size()-2;i++)
      {
        //d2Y[i]=1.0/(fabs(dY[i+1]-dY[i]))*(2.0*dY[i+1]-dY[i+2]/2.0-3.0/2.0*dY[i]);
        // формула f(x-h)-2f(x)+f(x+h)/h^2
        d2Y[i]=(y[i]-2*y[i+1]+y[i+2])/h/h;
        //d2Y[i]=(dY[i+2]-dY[i])/2.0/h;
      }
  return d2Y;
}

//-------------------------------------------------------------------
template <class T>
std::vector<T> calculateFirstDerivative(std::vector<T> & y, T h)
{
  std::vector<T> dY(y.size()-2);
  // Посчитаем производную методом конечных разностей
  // формула df/dx=(f(x+h)-f(x-h))/2/h;
  for(int i =0;i<y.size()-2;i++)
      {
        //dY[i]=1.0/(fabs(y[i+1]-y[i]))*(2.0*y[i+1]-y[i+2]/2.0-3.0/2.0*y[i]);
        dY[i]=(y[i+2]-y[i])/2.0/h;
      }
  return dY;
}





//---------------------------------------------------------------------------

template <class T>
void Rounding(T *pos, T* endPos)
{
    int S=1;
    for(;pos!=endPos;++pos)
    {
        int n=(int)(*pos*S)%10;
        if(n<5)
            *pos=floorl(*pos*S)/S;
        else
            *pos=ceill(*pos*S)/S;
    }
}

//--------------------------------------------------------------------------
// расчет среднего квадратичного отклонения
//--------------------------------------------------------------------------
template <class T>
T Sko (std::vector<T> const &x0,std::vector<T> const &x)
{
    int l=x0.size();
    T z=0;
    for(int i=0;i<l;i++)
    z+= pow(fabs((fabs(x[i])-fabs(x0[i]))),2);
    z/=(T)l;
    return sqrt(z);
}

//---------------------------------------------------
// расчет математического ожидания
//-------------------------------------------------

template <class T>
 T Mo (std::vector<T> const &x)
 {
    int l=x.size();
     T M=0;
     for(int i=0;i<l;i++)
     M+=x[i];
     return M/l;
 }

#endif
