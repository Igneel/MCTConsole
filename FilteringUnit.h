#ifndef FilteringUnitH
#define FilteringUnitH

#include <iterator>
#include "commonFunctions.h"

class FilterLowBand
{
public:
    FilterLowBand(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh);
    ~FilterLowBand();

    void FilterDataWithAutoShift(TSignal & inB,
    TSignal & inY,TSignal & outB,
    TSignal & outY);

    double FilterData (const std::vector<long double> &in,
    std::vector<long double> & out);

    void setFilterParams(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh);
    
private:

    void calculateImpulseResponse(unsigned int length,long double Fdisk, long double Fpropysk, long double Fzatyh);

    unsigned int N; //Длина фильтра и ведь она же должна считаться из частот)))
    long double Fd; //Частота дискретизации входных данных
    long double Fs; //Частота конца полосы пропускания 
    long double Fx; //Частота начала полосы затухания

    std::vector<long double> H;  //Импульсная характеристика фильтра
    std::vector<long double> H_id; //Идеальная импульсная характеристика
    std::vector<long double> W;   //Весовая функция
};

MyDataType BlockLowBandFilter(TSignal & inB,
TSignal & inY,TSignal & outB,
TSignal & outY,
size_t lengthFilter,MyDataType Fdisk,MyDataType Fpropysk,MyDataType Fzatyh,
size_t blockSize);

MyDataType TrForMassiveFilter(TSignal & inB,
TSignal & inY,TSignal & outB,
TSignal & outY,
int lengthFilter,MyDataType Fdisk,MyDataType Fpropysk,MyDataType Fzatyh);

double Filter (const std::vector<long double> &in,
std::vector<long double> & out, int length,
double Fdisk, double Fpropysk,double Fzatyh);

//---------------------------------------------------------------------------
#endif
