#ifndef INPUTFILEREADER_H
#define INPUTFILEREADER_H

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

/**
   \brief Class to encapsulate the parsing of input.txt files.
   
   As of 7/2/13 the inputs.txt file is arranged in the following manner:
   
   path to demands file
   number of weather files
   weather file 1
   weather file 2
   ...
   weather file n
   cell width
   quadrant angle
   deviation threshold
   node edge threshold
   angle offset
   lane width
   number of fixed nodes
   operational flexibility 1
   operational flexibility 2
   ...
   operational flexibility n
   timestamp

   When that format changes (like you decide that input.txt should be an XML file)
   encapsulate all of those changes here.
*/
class InputFileReader
{
private:
	std::string demandFile; /**< Path to demand file (i.e. demand.nom) */
    std::vector<std::string> weatherFileNames;  /**< Paths to weather files */
	double cellWidth; /**< Please describe this variable */
	double quadrantAngle;/**< Please describe this variable */
	double deviationThreshold; /**< Please describe this variable */
	double nodeEdgeThreshold; /**< Please describe this variable */
	double angleOffset; /**< Please describe this variable */
	double laneWidth; /**< Please describe this variable */
	unsigned demandShift; /**< Please describe this variable */
	unsigned int demandDrop; /**< Please describe this variable */
	int numFixedNodes; /**< Please describe this variable */
	std::vector<double> operFlex; /**< Please describe this variable */
	std::string timestamp; /**< When the input file was created */
public:
	InputFileReader(std::string inputFileName);
	InputFileReader(){}
public:
	std::string getDemandFile() const { return demandFile; }
	unsigned int getNumWeatherFiles() const { return weatherFileNames.size(); }
	std::string getWeatherFile(unsigned int index) const { return weatherFileNames[index]; }
	double getCellWidth() const { return cellWidth; }
	double getQuadrantAngle() const { return quadrantAngle; }
	double getDeviationThreshold() const { return deviationThreshold;}
	double getNodeEdgeThreshold() const { return nodeEdgeThreshold; }
	double getAngleOffset() const { return angleOffset; }
	double getLaneWidth() const { return laneWidth; }
	int getNumFixedNodes() const { return numFixedNodes; }
	unsigned int getDemandShift() const { return demandShift; }
	unsigned int getDemandDrop() const { return demandDrop; }
	std::vector<double> getOperFlex() const { return operFlex; }
	std::string getTimestamp() const { return timestamp; }
};

#endif // already included