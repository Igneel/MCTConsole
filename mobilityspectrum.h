#ifndef MOBILITYSPECTRUMH
#define MOBILITYSPECTRUMH

#include "commonFunctions.h"

using namespace std;

class mobilitySpectrum
{
private:
    int PointPerInt;

    typedef vector <long double> Data_spektr;

    typedef vector <long double> ImageDat;
    typedef vector <vector <long double> > mat;
    typedef vector <long double> Dat1;
    typedef vector <vector <long double> > Dat2;

    enum PeakType { MostSignificant, Significant, Slowdown};

    typedef vector< pair<long double, long double> > TLineSeries;

    unsigned long NumberOfPoints;
    int Power_spektr,GridPoints;

    Data_spektr MagField_spektr,GxxExp,GxyExp;

    ImageDat IntGxx;
    ImageDat IntGxy;
    ImageDat IntMagField;
    ImageDat Spectr_e;
    ImageDat Spectr_p;
    ImageDat Mobility;

    unsigned long SizeData;
    Dat1 B_spektr,Gxx_sp,Gxx_MC, Gxy_MC, Gxy_sp,Xr,Lv,Xv,Mv,Vpr;
    Dat2 Am,Qm,Cl,Cr,Cl_t,Cr_t,Cm,Cm_t          ;
    bool bulua;
    int  MSRight,MSLeft;

    long double W,F_s,A1,An,B1,Bn, Ves1, Ves2,Mu_max, Min_Spectr,Coef1,Coef2,Mu_min;

    TLineSeries electronMobilitySpectrum, holeMobilitySpectrum;

    TSignal resultElectronConductivity;
    TSignal resultHoleConductivity;
    TSignal resultMobility;

    std::vector<size_t> extremumHoleIndex;
    std::vector<size_t> extremumElectronIndex;

    std::vector<PeakType> holePeakType;
    std::vector<PeakType> electronPeakType;

    unsigned int  MaxPoints;

    void  InitArray();
    void  InitArray2();

    void  GetCoef(const Data_spektr &A,const Data_spektr &X, const long double b,
        long double &p0, long double &p1, long double &p2);
    void  GetLnLimits(int &Xmin, int &Xmax );

    void SetLength (vector<vector<long double> > &v,const size_t size1,const size_t size2);

    void  gram(const int N,const int m,const int l, Data_spektr & x, Data_spektr & f, mat & a);
    void  gauss(const int N, mat & a, Data_spektr & x);

    void  fi(const int n,const int m,const int l, Data_spektr & c, Data_spektr & x,
        const long double x1, long double &s);

    void  BAS(const int n,const int m,const int L,const long double x1, Data_spektr & x, Data_spektr & t);

    void  AddExpPoints(TLineSeries &ExpXX, TLineSeries &ExpXY);

    void  CS(Data_spektr& X,Data_spektr& F,Data_spektr& C,const long double p1,const long double pn);
    long double Sp(Data_spektr & X,Data_spektr &F, Data_spektr & C,const long double x1);
    void  MakeInterpolate(TLineSeries &Gxx, TLineSeries &Gxy,
        TLineSeries &ExpXX, TLineSeries &ExpXY);

    void  MakeMNK(const bool a, TLineSeries &Gxx, TLineSeries &Gxy, TLineSeries &ExpXX, TLineSeries &ExpXY);
    void  MakeLagranj();
    void  Tred2(const int n, Dat1 & d, Dat1 & e,
                     Dat2 & a, Dat2 & z, bool & fail);
    void  Imtql2(const int n,const long double macheps, Dat1 & d, Dat1 & e,
                     Dat2 & z, bool & fail);
    long double GetElem(const unsigned int j1, const int k1, const int i1);
    void  MakeMatrC();
    void  MakeMatrA();
    void  InverseMatrC(Dat2 & Ci,Dat2 & C,long double & Su,const int NP);
    long double S_s(const long double Mi);

    std::vector<long double> eigenValues;
    std::vector< std::vector<long double> > eigenVectors;


    size_t searchPeakRigthBorder(std::vector<long double> dh,std::vector<long double> d2h, size_t index);
    size_t searchPeakLeftBorder(std::vector<long double> dh,std::vector<long double> d2h, size_t index);
    void constructPeakCriteria(PeaksCriteria & peaksCriteria, TStringList * tsl, const std::vector<long double> & resMob, const std::vector<long double> & resCond, int index, int i, int j);

   void calculatePeakWeigth(PeaksCriteria & peaksCriteria,TStringList * tsl, 
  TSignal & resultMobility, TSignal & resultConductivity,
  std::vector<long double> & d, std::vector<long double> & d2, size_t extremumIndex);
    std::vector<long double> de;
    std::vector<long double> dh;
    std::vector<long double> d2e;
    std::vector<long double> d2h;

    std::vector <PeaksCriteria> vPC;
    AdditionalData additionalData;

public:

    std::vector <PeaksCriteria> getPeaksCriteria();
    AdditionalData getAdditionalData();

    size_t getResultSize();

    void saveEigenValues(std::string filename);
    void saveResults(std::string filename);
    void savePeakWeigth(std::string filename);

    long double calculatePeaksWeigth();
    long double calculatePeaksWeigth(std::string filename);
    void  MobilitySpectrumFunc(TLineSeries &LineSeries1, TLineSeries &Series5);

    long double getResultEY(const int i);

    long double getResultEX(const int i);

    long double getResultHY(const int i);

    long double getResultHX(const int i);

    std::vector<long double> getEigenValues();
    std::vector< std::vector<long double> > getEigenVectors();

    mobilitySpectrum(Data_spektr &MagneticFieldP, Data_spektr &Exx,
                     Data_spektr &Exy,const int size);

    void getExtremums(TSignal & holeConcentration, TSignal & holeMobility, TSignal & electronConcentration, TSignal & electronMobility);
    void getExtremumHoleIndex(std::vector<size_t> &v);
    void getExtremumElectronIndex(std::vector<size_t> &v);
    void getResults(TSignal & mobility, TSignal & HoleConductivity, TSignal & ElectronConductivity);

    long double calcConcentrationFromGp(long double G_p, long double Mu);
    size_t searchSignalSlowdown(TSignal & y, size_t startPosition, long double h);
    size_t searchSignificantPeak(TSignal & y, size_t startPosition, long double h);

};

#endif // MOBILITYSPECTRUM_H
