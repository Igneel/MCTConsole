//---------------------------------------------------------------------------

#ifndef NoiseUnitH
#define NoiseUnitH

#include <vector>
#include <cmath>
long double Mo (std::vector<long double> const &x, unsigned long l);
long double Sko (const long double *x0,const long double *x,unsigned long l);
void ShumAdding(std::vector<long double> const &x, std::vector<long double> &out,
                std::vector<long double> &ret, long double koeff, const unsigned long l);
void QuantumShumAdding(std::vector<long double> const &x, std::vector<long double> &out,
    std::vector<long double> &ret, const unsigned long l, long double mzr, int amplitude);

//---------------------------------------------------------------------------
#endif
