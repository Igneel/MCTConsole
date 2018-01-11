#ifndef DATASAVERH
#define DATASAVERH

#include "commonFunctions.h"



class DataSaver
{
public:
	DataSaver(std::string Temperature, std::string Current, std::string SampleInventoryNumber, std::string length, std::string width, std::string Thickness);
	~DataSaver();

	void setSampleDescription(std::string Temperature, std::string Current, std::string SampleInventoryNumber, std::string length, std::string width, std::string Thickness);
	void getSampleDescription(std::string &Temperature, std::string &Current, std::string &SampleInventoryNumber, std::string &length, std::string &width, std::string &Thickness);

	void SaveSampleDescription(std::string FileName);

void SaveData(DataKind dataKind,const TSignal *B,
const TSignal * HallEffect, const TSignal * MagnetoResistance,
SaveType saveType,std::string FileName);

    void setBaseFileName(std::string FileName);
	void setRoundNeeded(bool needRound);

	void setParamsType(ParamsType pt);

private:

	ParamsType paramsType;
	TSignal leftBound;
	TSignal rightBound;

	std::string BaseFileName;

	std::string defaultExtension;

	// Сохранение результатов и т.д.----------------------------------
	unsigned int NumberOfDecimalPlaces;
	bool isRoundNeeded;

	std::string Temperature; // температура образца во время измерений.
	std::string Current; // ток на образце, в амперах.

	std::string SampleInventoryNumber; // инвентарный номер образца.

	std::string SampleLength; // длина образца, мм
	std::string SampleWidth; // ширина образца, мм
	std::string SampleThickness; // толщина образца, мкм
    
    void SaveDataHelper(const TSignal *saveB,
	const TSignal * saveHall,
	const TSignal * saveResistance,SaveType mode,
	std::string FileName);
	inline void ReplaceDotsToComma(std::string &in, std::string & out);
	inline void ReplaceCommaToDots(std::string &in, std::string & out);
	
};

#endif
