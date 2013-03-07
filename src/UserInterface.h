#include "Quadrant.h"
#include "WeatherData.h"
#include "DemandProfile.h"
#include "RoutingDAG.h"
#include <string>

#define QUADRANT_NOT_GENERATED 0			// define if the quadrant is generated or not
#define QUADRANT_GENERATED 1

#define WEATHER_NOT_READ_IN 0				// define if the weather data files are read in
#define WEATHER_READ_IN 1

#define DEMAND_NOT_READ_IN 0				// define if the demand profile file is read in
#define DEMAND_READ_IN 1

#define ROUTINGDAG_NOT_GENERATED 0			// define if the routing DAG is generated
#define ROUTINGDAG_GENERATED 1

#define OPER_FLEX_GENERATED 0				// define if the operational flexibility is generated
#define OPER_FLEX_NOT_GENERATED 1

#define PI 3.14159265						// the value of PI

using namespace std;

#ifndef USERINTERFACE_H
#define USERINTERFACE_H

class UserInterface
{
public:
	UserInterface();						// constructor & destructor
	~UserInterface();
	/****************************************************************************************************************************/
	// general purpose functions
public:
	void reset();							// reset the status to a brand new routing instance
	void ProgramBegins();					// the project starts excuting from this function
private:
	void resetHelper();
	/****************************************************************************************************************************/
	// get user/file input or save to files functions
 private:
	string allInputs[200];
	int currentInput;  
	bool editQuadrant();					// ask for user input and edit the quadrant information
	bool readWeatherData();					// read in from weather data files 
	bool readDemandProfile();				// read in the demand profile information
	void saveTreeInformation();				// after generating the tree, save the information into an .xml file
	void inputOperationalFlexibility();		// ask user to input parameters related to operational flexity
	bool inputDemand();						// ask the user to input demand information
	void printQuadrantAndDemandInfo();		// use standard output to print the quadrant information on the screen
	bool inputDemandValid(string &input, int numDemands);		// tell if the user input demand is valid or not
	/****************************************************************************************************************************/
	// tree generating related functions
	bool generateTree();					// generate a bottommost robust tree
	bool tautenTree();						// tauten the generated bottommost tree
	/****************************************************************************************************************************/
private:
	Quadrant* quadrant;						// the quadrant in the routing context
	vector<WeatherData*> weatherDatas;		// a set of weather data files
	int numWeatherData;						// the number of weather data, equal to the size of the vector weatherDatas
	DemandProfile* demandProfile;			// demand profile information
	RoutingDAG* routingDAG;					// key component, the routing DAG structure
	vector<float> demandRNPs;				// the vector storing RNP values for each entry node
	float deviationThreshold;				// the threshold that a weather cell is considered hazardous
	float nodeEdgeThreshold;				// the threshold that we consider a thick edge or disk to be clear of weather
	/****************************************************************************************************************************/
	// control variables
private:
	int ctrl_QuadGenerated;					// mark if a quadrant was generated already
	int ctrl_WeatherReadIn;					// mark if a weather data file is already read in
	int ctrl_DemandReadIn;					// if the demand was read in already
	int ctrl_RoutingDAGGenerated;			// if a routingDAG is generated
	int ctrl_OperFlexGenerated;				// if the operational flxibility pairs are generated for an edge or a node
	/****************************************************************************************************************************/
	// variables related to coordinate system conversion between OPENGL(used in the windows version for displaying) and lati/long
private: 
	double centerLati,centerLong;			// the lati/long value of the center point (0, 0)
	double latiPerPixel, longPerPixel;		// lati/long that a single pixel represents
	string startTime, endTime;				// the time range of simulation
};

#endif
