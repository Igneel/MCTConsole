
#include "FilterParams.h"


FilterParams::FilterParams()
{
    SamplingFrequecy=10000;
    BandwidthFrequency=15;
    AttenuationFrequency=25;
    filterLength=50;
}

FilterParams::FilterParams(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,
                           MyDataType attenuationFrequency, long length)
{
    setFilterParams(samplingFrequecy, bandwidthFrequency, attenuationFrequency, length);
}

void FilterParams::setFilterParams(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,
                                   MyDataType attenuationFrequency, long length)
{
    SamplingFrequecy=samplingFrequecy;
    BandwidthFrequency=bandwidthFrequency;
    AttenuationFrequency=attenuationFrequency;
    filterLength=length;    
}

std::__cxx11::string FilterParams::getFilterParams()
{
    return "L"+std::to_string(filterLength)+"SF"+std::to_string(SamplingFrequecy)+"BF"+
    std::to_string(BandwidthFrequency)+"AF"+std::to_string(AttenuationFrequency);
}
