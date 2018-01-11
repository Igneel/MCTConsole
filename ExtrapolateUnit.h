#ifndef ExtrapolateUnitH
#define ExtrapolateUnitH

#include "commonFunctions.h"

bool testExtrapolateUnit();

long double determinant (long double ** Arr,int size);
long double determinant (std::vector<std::vector<long double> > & Arr,int size);

int curveFittingUniversal(long double * inX, long double *inY, const int length,
	long double *outKoef,const int powPolinom);
int curveFittingUniversal(std::vector<long double> * inX, std::vector<long double> *inY, 
	std::vector<long double> *outKoef,const int powPolinom);
void calculatePolinomByKoef(std::vector<long double> & inX, std::vector<long double> & koef, 
	std::vector<long double> & outF );

void LagrangeExtrapolation(TSignal & inX, TSignal & inY, TSignal & newX, TSignal & outY);

//---------------------------------------------------------------------------
#endif
