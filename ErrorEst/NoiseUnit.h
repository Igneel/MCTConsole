//---------------------------------------------------------------------------

#ifndef NoiseUnitH
#define NoiseUnitH

#include <vector>
#include <cmath>
//long double Mo (const long double *x,int l);
long double Mo (std::vector<long double> const &x,int l);
long double Sko (const long double *x0,const long double *x,int l);
//long double Sko (std::vector<long double> const &x0,std::vector<long double> const &x,int l);
//void ShumAdding(std::vector<long double> const *x,long double *out,long double *ret, long double koeff,const int l);
void ShumAdding(std::vector<long double> const &x,std::vector<long double> &out,std::vector<long double> &ret, long double koeff,const int l);

void QuantumShumAdding(std::vector<long double> const &x,std::vector<long double> &out,
    std::vector<long double> &ret,const int l,long double mzr, int amplitude);

//---------------------------------------------------------------------------
#endif
