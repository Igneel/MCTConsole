#ifndef FILTERPARAMSH
#define FILTERPARAMSH

#include "commonFunctions.h"

class FilterParams
{
public:
	FilterParams();
	FilterParams(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,MyDataType attenuationFrequency, int length);
    FilterParams(std::string samplingFrequecy,std::string bandwidthFrequency,std::string attenuationFrequency, std::string length);
	void setFilterParams(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,MyDataType attenuationFrequency, int length);

    std::string getFilterParams();

	MyDataType SamplingFrequecy;
	MyDataType BandwidthFrequency;
	MyDataType AttenuationFrequency;
	int filterLength;
};

#endif
