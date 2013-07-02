#ifndef INPUTFILEREADER_H
#define INPUTFILEREADER_H

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

class InputFileReader
{
private:
	std::vector<std::string> dataFileNames;
	std::string demandFile;
	double cellWidth;
	double quadrantAngle;
	double deviationThreshold;
	double nodeEdgeThreshold;
	double angleOffset;
	double laneWidth;
	unsigned int numFixedNodes;
	int operFlex1;
	int operFlex2;
	int operFlex3;
	std::string timestamp;
public:
	InputFileReader(std::string inputFileName);
	InputFileReader(){}
public:
	std::string getDemandFile() const { return demandFile; }
	unsigned int getNumWeatherFiles() const { return dataFileNames.size(); }
	std::string getWeatherFile(unsigned int index) const { return dataFileNames[index]; }
	double getCellWidth() const { return cellWidth; }
	double getQuadrantAngle() const { return quadrantAngle; }
	double getDeviationThreshold() const { return deviationThreshold;}
	double getNodeEdgeThreshold() const { return nodeEdgeThreshold; }
	double getAngleOffset() const { return angleOffset; }
	double getLaneWidth() const { return laneWidth; }
	unsigned int getNumFixedNodes() const { return numFixedNodes; }
	int getOperFlex1() const { return operFlex1; }
	int getOperFlex2() const { return operFlex2; }
	int getOperFlex3() const { return operFlex3; }
	std::string getTimestamp() const { return timestamp; }
};

#endif // already included