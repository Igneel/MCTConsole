//---------------------------------------------------------------------------
#include "FilteringUnit.h"

FilterLowBand::FilterLowBand(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh)
{
    N=0;
    calculateImpulseResponse(length, Fdisk, Fpropysk, Fzatyh);
}

FilterLowBand::~FilterLowBand()
{
    ;
}

void FilterLowBand::calculateImpulseResponse(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh)
{
    if (N!=length ||  Fd!=Fdisk || Fs!=Fpropysk || Fx != Fzatyh)
    {
        N = length; //Длина фильтра
        Fd = Fdisk; //Частота дискретизации входных данных 
        Fs = Fpropysk; //Частота конца полосы пропускания  
        Fx = Fzatyh; //Частота начала полосы затухания
         

        H.resize(N);
        H_id.resize(N);
        W.resize(N);
        //Расчет импульсной характеристики фильтра
        long double Fc = (Fs + Fx) / (2.0L * Fd);
        long double Wc = 2.0L*M_PIl*Fc;
        for (unsigned int i=0;i<N;++i)
        {
            if (i==N/2) H_id[i] = 2*Fc;
            else H_id[i] =2*Fc*sinl(Wc*(i-N/2))/(Wc*(i-N/2));
            // весовая функция Блекмена
            if (N>1)
            W[i] = 0.42L - 0.5L * cosl((2.0L*M_PIl*i) /( N-1.0L)) + 0.08L * cosl((4.0L*M_PIl*i) /( N-1.0L));
            else
            W[i]=1;
            H[i] = H_id[i] * W[i];
        }
        //Нормировка импульсной характеристики
        long double SUM=0.0L;
        for (unsigned int i=0; i<N; ++i) SUM +=H[i];
        for (unsigned int i=0; i<N; ++i) H[i]/=SUM; //сумма коэффициентов равна 1
    }
}

void FilterLowBand::setFilterParams(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh)
{
    calculateImpulseResponse(length, Fdisk, Fpropysk, Fzatyh);
}

double FilterLowBand::FilterData (const std::vector<long double> &in, std::vector<long double> & out)
{
    //Фильтрация входных данных
    unsigned int dataSize=in.size();
    out.resize(dataSize);
    for (unsigned int i=0; i<dataSize; ++i)
    {
        out[i]=0.0L;
        for (unsigned int j=0; j<=(i>N-1?N-1:i); ++j)// та самая формула фильтра
            out[i]+= H[j]*in[i-j];
    }
    return (N-1.0)/2.0;
}

void FilterLowBand::FilterDataWithAutoShift(TSignal & inB,
    TSignal & inY,TSignal & outB,
    TSignal & outY)
{
    auto lengthMassive=inY.size();
    if(lengthMassive==0)
    {
    return ;
    }
    outY.resize(lengthMassive);

    double k=FilterData(inY,outY); // вызываем фильтр
    k*=(max_elem(inB)-min_elem(inB))/(double)lengthMassive;// вычисляем сдвиг фаз
    // разность максимума и минимума на длину массива
    outB.resize(lengthMassive);
    for(auto i=0u;i<lengthMassive;i++) // выводим
    {
        outB[i]=inB[i]-k;
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// ВНИМАНИЕ!!! Напрямую не вызывать!!! Пользоваться трамплином!!!---------------
long double Filter (const std::vector<long double> &in, std::vector<long double> & out,
               unsigned int length, long double Fdisk, long double Fpropysk, long double Fzatyh)
{
unsigned int N = length; //Длина фильтра
long double Fd = Fdisk; //Частота дискретизации входных данных 2000
//long double Fs = Fpropysk; //Частота конца полосы пропускания  10
//long double Fx = Fzatyh; //Частота начала полосы затухания    20

std::vector<long double> H(N);  //Импульсная характеристика фильтра
std::vector<long double> H_id(N); //Идеальная импульсная характеристика
std::vector<long double> W(N);   //Весовая функция
    
//Расчет импульсной характеристики фильтра
long double Fc = (Fpropysk + Fzatyh) / (2.0L * Fd);
long double Wc = 2.0L * M_PIl * Fc;
/*
Строго говоря ширина перехода для функции Блэкмана фиксирована и составляет 5.5/N

Стало быть частота среза будет равна

*/
//Fzatyh= Fpropysk/Fd+5.5/(long double) N;
//Fc =  (Fpropysk + Fzatyh) / 2.0;
//Wc = 2*M_PI*Fc;

for (unsigned int i=0;i<N;++i)
{ // Идеальная импульсная характеристика фильтра нижних частот
if (i==N/2) H_id[i] = 2.0L * Fc;
else H_id[i] =2.0L *Fc*sinl(Wc*(i-N/2.0L))/(Wc*(i-N/2.0L));
}

for (unsigned int i=0;i<N;++i)
{
    // весовая функция Блекмена
    if(N>1)
        W [i] = 0.42L - 0.5L * cosl((2.0L*M_PIl*i) /(long double)N) + 0.08L * cosl((4.0L*M_PIl*i) /(long double)N);
    else
        W[i]=1.0L;
    H [i] = H_id[i] * W[i];
}

//Нормировка импульсной характеристики
long double SUM=0.0L;
for (unsigned int i=0; i<N; ++i) SUM +=H[i];
for (unsigned int i=0; i<N; ++i) H[i]/=SUM; //сумма коэффициентов равна 1


//Фильтрация входных данных
auto dataSize=in.size();
for (auto i=0u; i<dataSize; ++i)
{
out[i]=0.0L;
/*
for (unsigned int j=0; j<=(i>N-1?N-1:i); ++j)// та самая формула фильтра
out[i]+= H[j]*in[i-j];
}

*/
for (auto j=0u; j<N; ++j)// та самая формула фильтра
    if(i >= j)
        out[i]+= H[j]*in[i-j];
}
/*
TStringList * tsl= new TStringList;
for (int i=0; i<dataSize; ++i)
{
    tsl->Add(FloatToStr(out[i]));
}
tsl->SaveToFile("E:\\Дела\\Институт физики полупроводников\\Lcard\\MCTParamsExplorer\\filterOutPut.txt");
*/


return (N-1)/2.0;
}


long double TrForMassiveFilter(TSignal & in,TSignal &out,
int lengthFilter,long double Fdisk,
long double Fpropysk,long double Fzatyh)
{
    auto lengthMassive=in.size();
    if(lengthMassive==0)
    {
    return 0;
    }
    out.resize(lengthMassive);
    long double k=Filter(in,out,lengthFilter,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
    return k;
}


// функция трамплин для фильтра.  ----------------------------------------------
long double TrForMassiveFilter(std::vector<long double> & inB,
std::vector<long double> & inY,std::vector<long double> & outB,
std::vector<long double> & outY,int lengthFilter,long double Fdisk,
long double Fpropysk,long double Fzatyh)
{
    auto lengthMassive=inY.size();
    if(lengthMassive==0)
    {
    return 0;
    }
    outY.resize(lengthMassive);
    long double k=Filter(inY,outY,lengthFilter,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
    k*=(max_elem(inB)-min_elem(inB))/(long double)lengthMassive;// вычисляем сдвиг фаз
    // разность максимума и минимума на длину массива

//----------------------------------------------
//---------добавление для фильтрации магнитного поля
// ахтунг
//Filter(inB,outB,lengthFilter,Fdisk,Fpropysk,Fzatyh);
// внимание - она пока работает не корректно ибо фильтр наполняется только за
// его длину.
/*double *in2=new double[size];  // выделяем память
for(int i=0;i<size;i++)       // копируем
in2[i]=inS->XValues->Value[i];
double *out2=new double[size]; // выделяем память для выходных значений
double k2=Filter(in2,out2,size,length,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
k2*=(inS->XValues->MaxValue-inS->XValues->MinValue)/(double)inS->XValues->Count;// вычисляем сдвиг фаз
*/

//----------------------------------------------
    outB.resize(lengthMassive);
    for(auto i=0u;i<lengthMassive;i++) // выводим
    {
	    outB[i]=inB[i]-k;
    }
    return k;
}


MyDataType BlockLowBandFilter(TSignal & inB,
TSignal & inY,TSignal & outB,
TSignal & outY,
size_t lengthFilter,MyDataType Fdisk,MyDataType Fpropysk,MyDataType Fzatyh,
size_t blockSize)
{
if(blockSize<lengthFilter)
        lengthFilter=blockSize-1;

        TSignal tempInB;
        TSignal tempInS;
        TSignal tempOutB;
        TSignal tempOutS;
        size_t timesToRepeat=inB.size()/blockSize;
        for(unsigned int i=0; i <timesToRepeat;i++)
        {
                if(!outB.empty())
                {
                for(unsigned int j=0;j<lengthFilter;j++)
                {
                tempInB.push_back(outB.back());
                tempInS.push_back(outY.back());
                }
                }  

                for(unsigned int j=0+blockSize*i;j<(i+1)*blockSize;j++)
                {
                tempInB.push_back(inB[j]);
                tempInS.push_back(inY[j]);
                }
                tempOutB.resize(tempInB.size());
                tempOutS.resize(tempInB.size());
                TrForMassiveFilter(tempInB,tempInS,tempOutB,tempOutS,lengthFilter,
                Fdisk, Fpropysk, Fzatyh);
                for(unsigned int j=0;j<tempOutB.size();j++)
                {
                outB.push_back(tempOutB[j]);
                outY.push_back(tempOutS[j]);
                }
                tempOutB.clear();
                tempOutS.clear();
                tempInB.clear();
                tempInS.clear();
        }
        return 0;
}
