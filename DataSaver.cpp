#include "DataSaver.h"

//-------------------------------------------------------------------------------
DataSaver::DataSaver(std::string Temperature, std::string Current, std::string SampleInventoryNumber, 
	std::string length, std::string width, std::string Thickness)
{
	setSampleDescription(Temperature, Current, SampleInventoryNumber, length, width, Thickness);
	isRoundNeeded=true;
	defaultExtension=".txt";
	NumberOfDecimalPlaces=10;
	paramsType=DIRECT;

	leftBound.resize(3);
	leftBound[DIRECT]=0;
	leftBound[REVERSE]=-2;
	leftBound[COMBINE]=-2;
	rightBound.resize(3);
	rightBound[DIRECT]=2;
	rightBound[REVERSE]=0;
	rightBound[COMBINE]=2;
}
//-------------------------------------------------------------------------------
void DataSaver::setSampleDescription(std::string Temperature, std::string Current, 
	std::string SampleInventoryNumber, std::string length, std::string width, std::string Thickness)
{
	this->Current=Current;
	this->SampleInventoryNumber=SampleInventoryNumber;
	this->Temperature=Temperature;	
	this->SampleLength=length;
	this->SampleWidth=width;
	this->SampleThickness=Thickness;
}
//-------------------------------------------------------------------------------
void DataSaver::getSampleDescription(std::string &Temperature, std::string &Current, 
	std::string &SampleInventoryNumber, std::string &length, std::string &width, std::string &Thickness)
{

	Current=this->Current;
	SampleInventoryNumber=this->SampleInventoryNumber;
	Temperature=this->Temperature;
	length=this->SampleLength;
	width=this->SampleWidth;
	Thickness=this->SampleThickness;
}
//-------------------------------------------------------------------------------
void DataSaver::SaveSampleDescription(std::string FileName)
{
	std::string NewFileName;
	// Нам пока не нужны уникальные имена. Нам нужны фиксированные.
    //NewFileName=FileName.substr(0,FileName.length()-3)+"Description"+defaultExtension;
    NewFileName=FileName.substr(0,FileName.length()-3)+defaultExtension;
	TStringList * tsl=new TStringList();

    tsl->push_back("SampleInventoryNumber\t"+SampleInventoryNumber);
    tsl->push_back("Temperature\t"+Temperature);
    tsl->push_back("Current\t"+Current);
    tsl->push_back("SampleLength\t"+SampleLength);
    tsl->push_back("SampleWidth\t"+SampleWidth);
    tsl->push_back("SampleThickness\t"+SampleThickness);

	//FileName+=defaultExtension;

    SaveToFile(tsl,NewFileName);

	delete tsl;
}
//-------------------------------------------------------------------------------
DataSaver::~DataSaver()
{
	;
}

//-------------------------------------------------------------------------------

void DataSaver::SaveData(DataKind dataKind,const TSignal *B,
const TSignal * HallEffect, const TSignal * MagnetoResistance,
SaveType saveType,std::string FileName)
{

	std::string NewFileName;
	NewFileName=FileName+"SN"+SampleInventoryNumber+"T"+Temperature+"I"+Current;
	switch(dataKind)
	{
	    case CURRENT_DATA:
	        NewFileName+="CurrD";
	        break;
	    case FILTERED_DATA:
	        NewFileName+="FiltD";
	        break;
	    case EXTRAPOLATED_DATA:
	        NewFileName+="ExtrD";
	        break;
	    case AVERAGED_DATA:
	    	NewFileName+="AvD";
	    	break;
	    case TENZOR_DATA:
	    	NewFileName+="TenzorD";
	    	break;
	    case MOBILITY_DATA:
	    	NewFileName+="MobSpecD";
	    	break;
	    case FIT_DATA:
	    	NewFileName+="FitD";
	    	break;
	}

	SaveDataHelper(B,HallEffect,MagnetoResistance,saveType,NewFileName);

}
//-------------------------------------------------------------------------------

void DataSaver::SaveDataHelper(const TSignal *saveB,
const TSignal * saveHall,
const TSignal * saveResistance,SaveType mode, std::string FileName)
{
    TStringList * tsl=new TStringList();

    TSignal savingXData(saveB->begin(),saveB->end());
    TSignal savingY1Data(saveResistance->begin(),saveResistance->end());
    TSignal savingY2Data(saveHall->begin(),saveHall->end());

    auto length=savingXData.size();
    if(isRoundNeeded==true)
	{
        roundM(savingXData,NumberOfDecimalPlaces);
        roundM(savingY1Data,NumberOfDecimalPlaces);
        roundM(savingY2Data,NumberOfDecimalPlaces);
	}

	if (mode==POINTS_11) 
	{
        FileName+="11P";
        TSignal tempX;
        TSignal tempY;
        thiningSignal(savingXData,savingY1Data,tempX,tempY,leftBound[paramsType],rightBound[paramsType],11);
        savingY1Data = tempY;

        thiningSignal(savingXData,savingY2Data,tempX,tempY,leftBound[paramsType],rightBound[paramsType],11);
        savingXData = tempX;
        savingY2Data = tempY;

        for (unsigned int i = 0; i < savingXData.size(); ++i)
        {
            tsl->push_back(std::to_string(savingXData[i])+"\t"+
            std::to_string(savingY1Data[i])+"\t"+
            std::to_string(savingY2Data[i]));
        }
	}
    if (mode==POINTS_21)
    {
        FileName+="21P";
        
        TSignal tempX;
        TSignal tempY;
        thiningSignal(savingXData,savingY1Data,tempX,tempY,leftBound[paramsType],rightBound[paramsType],21);
        savingY1Data = tempY;

        thiningSignal(savingXData,savingY2Data,tempX,tempY,leftBound[paramsType],rightBound[paramsType],21);
        savingXData = tempX;
        savingY2Data = tempY;

        for (unsigned int i = 0; i < savingXData.size(); ++i)
        {
            tsl->push_back(std::to_string(savingXData[i])+"\t"+
            std::to_string(savingY1Data[i])+"\t"+
            std::to_string(savingY2Data[i]));
        }        
    }
	if(mode==ALL_POINTS)
	{
        FileName+="AllP";
        for(auto i=0u;i<length;i++)
		{
            tsl->push_back(std::to_string(savingXData[i])+"\t"+
            std::to_string(savingY1Data[i])+"\t"+
            std::to_string(savingY2Data[i]));
		}
	}
    //std::string text=tsl->Text.c_str();
   
    //ReplaceCommaToDots(text,text);
    //tsl->Text=text.c_str();

    FileName+=defaultExtension;
    if(tsl->size()>0){
        SaveToFile(tsl,FileName);
    }

	delete tsl;
}
//-------------------------------------------------------------------------------
void DataSaver::setRoundNeeded(bool needRound)
{
    isRoundNeeded=needRound;
}

//-------------------------------------------------------------------------------
inline void DataSaver::ReplaceDotsToComma(std::string &in, std::string & out)
{
    unsigned int findIndex=0;
	std::string s=in;
	std::string strToReplaceWhich=","; // на что меняем
	std::string strToSearch=".";   // что ищем
	while ((findIndex=s.find(strToSearch,findIndex))!=std::string::npos)
    {
        s.replace(s.begin()+findIndex,s.begin()+findIndex+strToSearch.length(),
        strToReplaceWhich.begin(),strToReplaceWhich.end());
    } 
    out=s;
}
//-------------------------------------------------------------------------------
inline void DataSaver::ReplaceCommaToDots(std::string &in, std::string & out)
{
    unsigned int findIndex=0;
	std::string s=in;
	std::string strToReplaceWhich="."; // на что меняем
	std::string strToSearch=",";   // что ищем
    while ((findIndex=s.find(strToSearch,findIndex))!=std::string::npos)
    {
    	s.replace(s.begin()+findIndex,s.begin()+findIndex+strToSearch.length(),
        strToReplaceWhich.begin(),strToReplaceWhich.end());
    } 
    out=s;
}

void DataSaver::setParamsType(ParamsType pt)
{
	paramsType=pt;
}
