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
   quadrant angular width
   lane width
   number of fixed nodes
   operational flexibility 1
   operational flexibility 2
   ...
   operational flexibility n
   outputFilename

   When that format changes (like you decide that input.txt should be an XML file)
   encapsulate all of those changes here.
*/
class InputFileReader
{
 private:
  std::string demandFile; /**< Path to demand file (i.e. demand.nom) */
  std::vector<std::string> weatherFileNames;  /**< Paths to weather files */
  double cellWidth; /**< Please describe this variable */
  double quadrantAngle;/**< Starting angle of our wedge. */
  /** Ignore weather that makes pilots less than this likely to fly around */
  double deviationThreshold; 
  double nodeEdgeThreshold; /**< 1-sum(prob emsembles impacted) < this */
  double angularWidth; /**< Angular width of quadrant */
  /** Clearance required from weather (in nm)
      Set non-zero to override positive demand RNPs from demand file. */
  double laneWidth; 
  unsigned int demandShift; /**< Please describe this variable */
  unsigned int demandDrop; /**< Please describe this variable */
  /** Maximum number of fixed nodes to generate
      Negative value used to indicate no limit */
  int numFixedNodes; 
  std::vector<double> operFlex; /**< Please describe this variable */
  std::string outputFilename; /**< Name of output XML file */
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
  double getAngularWidth() const { return angularWidth; }
  double getLaneWidth() const { return laneWidth; }
  int getNumFixedNodes() const { return (int)numFixedNodes; }
  unsigned int getDemandShift() const { return demandShift; }
  unsigned int getDemandDrop() const { return demandDrop; }
  const std::vector<double> &getOperFlex() const { return operFlex; }
  std::string getOutputName() const { return outputFilename; }
 private:
  bool areInputsLegal();
};

#endif // already included
