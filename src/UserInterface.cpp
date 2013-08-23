// #define LEARNING_ABOUT_TREE

#include "UserInterface.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

/*****************************************/
// header files from within the project
#include "NodeAndEdge.h"
#include "InputFileReader.h"

#if !defined(DO_NOT_CONVERT_LAT_LON_TO_PIXELS)
#  include "RoutingDAG.h"
#endif

#if defined(SUPPRESS_OUTPUT)
#define cout ostream(0).flush()
#endif


/* Have to update this if we ever allow more than 4 as a parameter for demand_drop
There is also a 4 that needs to be changed down below, in the declaration of an array called temp_demands[4]
There is a ruby file (combination.rb) to generate this line */  

int combinations[][4] = { {0,-1,-1,-1,}, {1,-1,-1,-1,}, {2,-1,-1,-1,}, {3,-1,-1,-1,}, {4,-1,-1,-1,}, {5,-1,-1,-1,}, {6,-1,-1,-1,}, {7,-1,-1,-1,}, {8,-1,-1,-1,}, {0,1,-1,-1,}, {0,2,-1,-1,}, {0,3,-1,-1,}, {0,4,-1,-1,}, {0,5,-1,-1,}, {0,6,-1,-1,}, {0,7,-1,-1,}, {0,8,-1,-1,}, {1,2,-1,-1,}, {1,3,-1,-1,}, {1,4,-1,-1,}, {1,5,-1,-1,}, {1,6,-1,-1,}, {1,7,-1,-1,}, {1,8,-1,-1,}, {2,3,-1,-1,}, {2,4,-1,-1,}, {2,5,-1,-1,}, {2,6,-1,-1,}, {2,7,-1,-1,}, {2,8,-1,-1,}, {3,4,-1,-1,}, {3,5,-1,-1,}, {3,6,-1,-1,}, {3,7,-1,-1,}, {3,8,-1,-1,}, {4,5,-1,-1,}, {4,6,-1,-1,}, {4,7,-1,-1,}, {4,8,-1,-1,}, {5,6,-1,-1,}, {5,7,-1,-1,}, {5,8,-1,-1,}, {6,7,-1,-1,}, {6,8,-1,-1,}, {7,8,-1,-1,}, {0,1,2,-1,}, {0,1,3,-1,}, {0,1,4,-1,}, {0,1,5,-1,}, {0,1,6,-1,}, {0,1,7,-1,}, {0,1,8,-1,}, {0,2,3,-1,}, {0,2,4,-1,}, {0,2,5,-1,}, {0,2,6,-1,}, {0,2,7,-1,}, {0,2,8,-1,}, {0,3,4,-1,}, {0,3,5,-1,}, {0,3,6,-1,}, {0,3,7,-1,}, {0,3,8,-1,}, {0,4,5,-1,}, {0,4,6,-1,}, {0,4,7,-1,}, {0,4,8,-1,}, {0,5,6,-1,}, {0,5,7,-1,}, {0,5,8,-1,}, {0,6,7,-1,}, {0,6,8,-1,}, {0,7,8,-1,}, {1,2,3,-1,}, {1,2,4,-1,}, {1,2,5,-1,}, {1,2,6,-1,}, {1,2,7,-1,}, {1,2,8,-1,}, {1,3,4,-1,}, {1,3,5,-1,}, {1,3,6,-1,}, {1,3,7,-1,}, {1,3,8,-1,}, {1,4,5,-1,}, {1,4,6,-1,}, {1,4,7,-1,}, {1,4,8,-1,}, {1,5,6,-1,}, {1,5,7,-1,}, {1,5,8,-1,}, {1,6,7,-1,}, {1,6,8,-1,}, {1,7,8,-1,}, {2,3,4,-1,}, {2,3,5,-1,}, {2,3,6,-1,}, {2,3,7,-1,}, {2,3,8,-1,}, {2,4,5,-1,}, {2,4,6,-1,}, {2,4,7,-1,}, {2,4,8,-1,}, {2,5,6,-1,}, {2,5,7,-1,}, {2,5,8,-1,}, {2,6,7,-1,}, {2,6,8,-1,}, {2,7,8,-1,}, {3,4,5,-1,}, {3,4,6,-1,}, {3,4,7,-1,}, {3,4,8,-1,}, {3,5,6,-1,}, {3,5,7,-1,}, {3,5,8,-1,}, {3,6,7,-1,}, {3,6,8,-1,}, {3,7,8,-1,}, {4,5,6,-1,}, {4,5,7,-1,}, {4,5,8,-1,}, {4,6,7,-1,}, {4,6,8,-1,}, {4,7,8,-1,}, {5,6,7,-1,}, {5,6,8,-1,}, {5,7,8,-1,}, {6,7,8,-1,}, {0,1,2,3,}, {0,1,2,4,}, {0,1,2,5,}, {0,1,2,6,}, {0,1,2,7,}, {0,1,2,8,}, {0,1,3,4,}, {0,1,3,5,}, {0,1,3,6,}, {0,1,3,7,}, {0,1,3,8,}, {0,1,4,5,}, {0,1,4,6,}, {0,1,4,7,}, {0,1,4,8,}, {0,1,5,6,}, {0,1,5,7,}, {0,1,5,8,}, {0,1,6,7,}, {0,1,6,8,}, {0,1,7,8,}, {0,2,3,4,}, {0,2,3,5,}, {0,2,3,6,}, {0,2,3,7,}, {0,2,3,8,}, {0,2,4,5,}, {0,2,4,6,}, {0,2,4,7,}, {0,2,4,8,}, {0,2,5,6,}, {0,2,5,7,}, {0,2,5,8,}, {0,2,6,7,}, {0,2,6,8,}, {0,2,7,8,}, {0,3,4,5,}, {0,3,4,6,}, {0,3,4,7,}, {0,3,4,8,}, {0,3,5,6,}, {0,3,5,7,}, {0,3,5,8,}, {0,3,6,7,}, {0,3,6,8,}, {0,3,7,8,}, {0,4,5,6,}, {0,4,5,7,}, {0,4,5,8,}, {0,4,6,7,}, {0,4,6,8,}, {0,4,7,8,}, {0,5,6,7,}, {0,5,6,8,}, {0,5,7,8,}, {0,6,7,8,}, {1,2,3,4,}, {1,2,3,5,}, {1,2,3,6,}, {1,2,3,7,}, {1,2,3,8,}, {1,2,4,5,}, {1,2,4,6,}, {1,2,4,7,}, {1,2,4,8,}, {1,2,5,6,}, {1,2,5,7,}, {1,2,5,8,}, {1,2,6,7,}, {1,2,6,8,}, {1,2,7,8,}, {1,3,4,5,}, {1,3,4,6,}, {1,3,4,7,}, {1,3,4,8,}, {1,3,5,6,}, {1,3,5,7,}, {1,3,5,8,}, {1,3,6,7,}, {1,3,6,8,}, {1,3,7,8,}, {1,4,5,6,}, {1,4,5,7,}, {1,4,5,8,}, {1,4,6,7,}, {1,4,6,8,}, {1,4,7,8,}, {1,5,6,7,}, {1,5,6,8,}, {1,5,7,8,}, {1,6,7,8,}, {2,3,4,5,}, {2,3,4,6,}, {2,3,4,7,}, {2,3,4,8,}, {2,3,5,6,}, {2,3,5,7,}, {2,3,5,8,}, {2,3,6,7,}, {2,3,6,8,}, {2,3,7,8,}, {2,4,5,6,}, {2,4,5,7,}, {2,4,5,8,}, {2,4,6,7,}, {2,4,6,8,}, {2,4,7,8,}, {2,5,6,7,}, {2,5,6,8,}, {2,5,7,8,}, {2,6,7,8,}, {3,4,5,6,}, {3,4,5,7,}, {3,4,5,8,}, {3,4,6,7,}, {3,4,6,8,}, {3,4,7,8,}, {3,5,6,7,}, {3,5,6,8,}, {3,5,7,8,}, {3,6,7,8,}, {4,5,6,7,}, {4,5,6,8,}, {4,5,7,8,}, {4,6,7,8,}, {5,6,7,8,} };

// the constructor sets all the control variables to their initial status
UserInterface::UserInterface()
{
	// initialize ctrl variables to default and generate necassary objects
	resetHelper();	
}

// the destructor reclaims all the allocated memory to avoid memory leaks
UserInterface::~UserInterface()
{
	if(ctrl_QuadGenerated == QUADRANT_GENERATED)			delete quadrant;
	if(ctrl_WeatherReadIn == WEATHER_READ_IN)				delete demandProfile;
	if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)	delete routingDAG;
}

// used to generate a brand new tree generating instance, everything is reset to original value
void UserInterface::reset()
{
	// deleting the existing objects
	if(ctrl_QuadGenerated == QUADRANT_GENERATED)
	{
		delete quadrant;
	}
	if(ctrl_WeatherReadIn == WEATHER_READ_IN)
	{ 
		delete demandProfile;
	}
	if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)	
	{
		delete routingDAG;
	}
	weatherDataSets.clear();
	demandRNPs.clear();
	resetHelper(); // set all the ctrl variables to default and generate new objects
}

// set the system to the starting status, including variables and object initializations
void UserInterface::resetHelper()
{
	ctrl_QuadGenerated = QUADRANT_GENERATED;				// there is always a default quadrant
	ctrl_WeatherReadIn = WEATHER_NOT_READ_IN;				// current default: no weather data is read in
	ctrl_DemandReadIn = DEMAND_NOT_READ_IN;					// and there is no demand profile read in
	ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;	// the routing dag is not generated at the beginning
	ctrl_OperFlexGenerated = OPER_FLEX_NOT_GENERATED;		// the operational flxitibility pairs are not generated at the beginning

	// generate the set of necessary objects
	quadrant = new Quadrant();
	demandProfile = new DemandProfile();
	routingDAG = new RoutingDAG();
	deviationThreshold = (double) 0.7;
	nodeEdgeThreshold = (double) 0.8; // default values of the thresholds
}


/****************************************************************************************/
// RTPrototype Functions

// The project starts excuting here, called by the main function
bool UserInterface::ProgramBegins(std::string inputFile)
{

	inputs = InputFileReader(inputFile);
	outputFileName = inputs.getOutputName();

	char userInput = '0';	     // the user input will always be a number, deciding what to do next
	bool startOver = false;    // when the user chooses to do everything again
	std::cout << "Robust Tree Generator:\n";
	std::cout << "The software is used to generate robust trees given the demand profile and weather data information." << std::endl;

	if(startOver) // when doing everything again, reset everything first
	{
		reset();
		startOver = false;
	}
	if(!readDemandProfile()) // at the start, read in demand profile
	{
		std::cout << "\nFailed to read in the demand profile, double check the demand path listed." << std::endl;
		return false;
	}
	std::cout << "Demand profile successfully read in." << std::endl;
	if(!readWeatherData())
	{
		std::cout << "\nFailed to successfully read in the weather data, double check the weather directories." << std::endl;
		return false;
	}
	std::cout << "\nWeather files are succesfully read in!" << std::endl;
	(*quadrant).setAngle( inputs.getQuadrantAngle() );
	(*quadrant).setAngularWidth( inputs.getQuadrantAngle() );
	std::cout << "Current angle: " << quadrant->getAngle() << std::endl;

	return makeDAG();
}

void UserInterface::makeRTPTreeAndFinish()
{
	std::cout << "Generating tree." << std::endl;
	generateTree();

	std::cout << "Tautening tree." << std::endl;
	tautenTree();

	std::cout << "\nDoing operational Flexibility stuff." << std::endl;
	inputOperationalFlexibility();
	std::cout << "\nSaving tree information." << std::endl;

	saveTreeInformation();

}

// print the information of the quadrant and the demand rnps for each entry node
void UserInterface::printQuadrantAndDemandInfo()
{
	if(ctrl_QuadGenerated == QUADRANT_GENERATED)  
	{
		// first print the current quadrant information and demand information
		std::cout << "\nThe center of the quadrant is currently at lati/long ("<<setprecision(5)<<centerLati+quadrant->getCenterX()*latiPerPixel<<", ";
		std::cout << setprecision(5)<<centerLong+quadrant->getCenterY()*longPerPixel<<").";
		std::cout << "\nThe inner radius (in nm) is "<<setprecision(5)<<quadrant->getiRadius()*NMILESPERPIXEL<<", the outer radius is ";
		std::cout << setprecision(5)<<quadrant->getoRadius()*NMILESPERPIXEL<<".";
		std::cout << "\nThe inner altitude (in feet) is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getiHeight()*ALTITUDE_PER_PIXEL<<", ";
		std::cout << "the outer altitude is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getoHeight()*ALTITUDE_PER_PIXEL<<".";
		std::cout << "\nThe angle of the right boundary of the quadrant is (relative to +x axis) "<<setprecision(5)<<quadrant->getAngle()*180/PI<<" degrees."<<std::endl;
		// then print the demand information
		std::cout << "\nCurrently there are "<<demandRNPs.size()<<" entry points evenly spaced on the outer boundary of the quadrant.";
		if(!demandRNPs.empty())
		{
			std::cout << "\nthe rnp requirements are (in nm, from right to left) ";
			for(unsigned i=0; i<demandRNPs.size(); i++)
				std::cout << setprecision(5)<<demandRNPs[i]*NMILESPERPIXEL<<"  ";
		}
	}
	else
	{
		std::cout << "\nThe quadrant is not generated yet!"<<std::endl;
	}
}

/*****************************************************************************************************************************************/

/** 
\brief Generates a tree
*/
bool UserInterface::makeDAG()
{
	double quadrantAngularWidth;
	double lane_width;
	int maxFixNodes;

	// only when the quadrant is generated and weather data is read in can we generate a tree
	if(ctrl_QuadGenerated == QUADRANT_GENERATED && ctrl_WeatherReadIn == WEATHER_READ_IN && ctrl_DemandReadIn == DEMAND_READ_IN)
	{
		double minDistBetweenMergeNodes = 5/NMILESPERPIXEL; // WILLXYZ make this read from config file

		std::cout << "Reading in deviation threshold." << std::endl;
		deviationThreshold = inputs.getDeviationThreshold();

		std::cout << "Reading in nodeEdgeThreshold." << std::endl;
		nodeEdgeThreshold = inputs.getNodeEdgeThreshold();
		routingDAG->setminimumDistanceBetweenMergingNodes(minDistBetweenMergeNodes);
		std::cout <<  std::endl << "Generating a bottommost merge tree, please wait..." << std::endl;
		/************************************************************************************************/
		ctrl_OperFlexGenerated = OPER_FLEX_NOT_GENERATED;				// when generating a new tree, the Oper-Flex pairs need to be generated again
		// first, generate the entry and fixed nodes, then the internal nodes
		routingDAG->reset();											// a brand new routing instance
		std::cout << std::endl << "Finished resetting edges." << std::endl;

		quadrantAngularWidth = inputs.getAngularWidth();
		lane_width = inputs.getLaneWidth();
		// Note: if lane width is negative, we will use lane widths from 
		// demand files.
		if (lane_width > 0) 
		{
			for (unsigned int i = 0; i < demandRNPs.size(); i++) 
			{
				if (demandRNPs[i] > 0) 
				{
					demandRNPs[i] = lane_width;
				}
			}
		}
#if !defined(DO_NOT_CONVERT_LAT_LON_TO_PIXELS)
		std::cout << "Demand RNPS (in pixels): ";
#else
		std::cout << "Demand RNPS (in nm): ";
#endif
		for(unsigned int i = 0; i < demandRNPs.size(); i++) 
		{
			std::cout << demandRNPs[i] << " ";
		}
		std::cout << std::endl;

		std::cout << "Max fixed nodes: " << inputs.getNumFixedNodes() << std::endl;
		maxFixNodes = inputs.getNumFixedNodes(); // a negative parameter is used sometimes to indicate no limit
		std::cout << "START GENERATING DAG" << std::endl;
	}
	else // the weather data and demand profile have to be read in first
	{
		std::cerr << "\nPlease read in or generate the demand profile and weather data first." << std::endl;
		return false;
	}

	if(quadrant->generateDAG(demandRNPs, demandRNPs.size(), deviationThreshold, nodeEdgeThreshold, weatherDataSets, routingDAG, quadrantAngularWidth, maxFixNodes))
	{
		std::cout << "FINISHED DAG" << std::endl;
		std::cout << std::endl << "Generating edge set..." << std::endl;
		routingDAG->generateEdgeSet();								// generate the edges in the searching DAG
		std::cout << std::endl << "Edges (Routing DAG) generated." << std::endl;
		ctrl_RoutingDAGGenerated = ROUTINGDAG_GENERATED;
		std::cout << std::endl << "Generating Tree..." << std::endl;
		// generate the tree here
		// =======
		std::cout << "START GENERATING TREE" << std::endl << "Demand RNPs: ";
		for (unsigned int i = 0; i < demandRNPs.size(); i++) 
		{
			std::cout << demandRNPs[i] << "  ";
		}
		std::cout << std::endl;
	}
	else 
	{ 
		std::cout << "Failed to generate the DAG."; 
		return false;
	}
	return true;
}
	bool UserInterface::generateTree()
	{

	unsigned int demand_shift = inputs.getDemandShift();
	unsigned int demand_drop  = inputs.getDemandDrop();

	if(!routingDAG->generateTree(weatherDataSets, demandRNPs, deviationThreshold, nodeEdgeThreshold))
	{
		// Will added some code 04/2013 to add a demand shifting scheme to hopefully reduce tree generation failure
		bool demand_shift_success = false;
		if (demand_shift == 1) 
		{
			std::cout << std::endl << "Beginning demand shifting." << std::endl;
			std::cout << "Variables: demand_shift " << demand_shift << "; demand drop " << demand_drop << "; demand_shift_success " << demand_shift_success << std::endl;
			/* BEWARE: POTENTIAL BUG (due to Will not being certain what Shang's code does)
			It is possible that you can only attempt to generate a tree once before needing to redo everything.
			If this is the case, we're in trouble.

			It shouldn't be an issue, there is a 'resetTree()' method called at the beginning of the definition of 
			generateTree(), but it is a possible concern.
			*/

			/* We will try two techniques to make the tree more robust:
			(1) We will try shifting the demand up one or down one, and then generating the tree. If that fails, we will try
			(2) Deleting one demand node at a time, and then try generating the tree.

			The next while loop is part (1)*/
			for (unsigned int demand_shift_index = 0; 
				demand_shift_index < demandRNPs.size() && !demand_shift_success;
				demand_shift_index++) 
			{

				if (demand_shift_index > 0 && demandRNPs[demand_shift_index - 1] == 0 && !demand_shift_success) 
				{
					// Swap the values of demandRNPs at the current index and the one below
					demandRNPs[demand_shift_index - 1] = demandRNPs[demand_shift_index];
					demandRNPs[demand_shift_index] = 0;

					// If the tree is generated successfully, exit this block gracefully
					if (routingDAG->generateTree(weatherDataSets, demandRNPs, deviationThreshold, nodeEdgeThreshold)) 
					{
						demand_shift_success = true;
						std::cout << std::endl << "A downward swap generated a tree on the following index: " << demand_shift_index << std::endl;
					}

					// Swap them back
					demandRNPs[demand_shift_index] = demandRNPs[demand_shift_index - 1];
					demandRNPs[demand_shift_index - 1] = 0;

				}
				if (demand_shift_index + 1 < demandRNPs.size() && demandRNPs[demand_shift_index + 1] == 0 && !demand_shift_success) 
				{
					// Swap the values of demandRNPs at the current index and the one above
					demandRNPs[demand_shift_index + 1] = demandRNPs[demand_shift_index];
					demandRNPs[demand_shift_index] = 0;

					// If the tree is generated successfully, exit this block gracefully
					if (routingDAG->generateTree(weatherDataSets, demandRNPs, deviationThreshold, nodeEdgeThreshold)) {
						demand_shift_success = true;
						std::cout << std::endl << "An upward swap generated a tree on the following index: " << demand_shift_index << std::endl;
					}

					// Swap them back
					demandRNPs[demand_shift_index] = demandRNPs[demand_shift_index + 1];
					demandRNPs[demand_shift_index + 1] = 0;
				}
			}
		}
		if (demand_drop >= 1 && !demand_shift_success) 
		{
			/* Begin phase (2) of the demand shifting, aka demand_drop */
			std::cout << std::endl <<  "Demand shifting unsuccessful, beginning demand dropping." << std::endl;

			std::vector<double> temp_demands(4); // This 4 needs to be changed if combinations is ever changed -- HARDCODED
			// the below 255 is the length of ``combinations`` -- HARDCODED
			for (unsigned int k = 0, comb_i = 0; 
				k < demand_drop && comb_i < 255 && !demand_shift_success;
				comb_i++ ) 
			{
				for (int j = 0; j < 4; ++j) 
				{
					temp_demands[j] = 0;
				}

				// std::cout << std::endl;
				for (int j = 0; j < 4; ++j) 
				{
					if(combinations[comb_i][j] <   (int) 0 ||
						combinations[comb_i][j] >= (int)demandRNPs.size())
					{
						continue;
					}
					if (demandRNPs[combinations[comb_i][j]] > 0) 
					{
						temp_demands[j] = demandRNPs[ combinations[comb_i][j] ] ;
						demandRNPs[combinations[comb_i][j]] = 0;
					}
					// HARDCODED VALUE
					if (j == 3) 
					{
						if (routingDAG->generateTree(weatherDataSets, demandRNPs, deviationThreshold, nodeEdgeThreshold)) 
						{
							demand_shift_success = true;
							std::cout << std::endl << "A demand drop generated a tree on the following indices: ";
							for(unsigned int i = 0; i <= k; i++) 
							{
								std::cout << combinations[comb_i][i] << " ";
							}
							std::cout << std::endl;

							std::cout << std::endl << "Demand RNPs after the drop:" << std::endl;
							for(unsigned int i = 0; i < demandRNPs.size(); i++) 
							{
								std::cout << demandRNPs[i] << " ";
							}
							std::cout << std::endl;
						}
						else 
						{
							for (int jj = 0; jj < 4; ++jj) 
							{
								demandRNPs[combinations[comb_i][jj]] = temp_demands[jj];
								if (demandRNPs[combinations[comb_i][jj]] < 0.00000001) 
								{
									demandRNPs[combinations[comb_i][jj]] = 0;
								}
							}
						}
					}
				}
				// Find the new k value, which is the last positive number in the array -- HARDCODED
				for (int j = 0; j < 4; ++j) 
				{
					if (combinations[comb_i][j] < 0) 
					{
						k = j - 1;
						break;
					}
					if (j == 3)
					{ // j would go on to being 4, but can't because of the loop condition
						k = 3; 
					}

				}
			}
		}
		if (!demand_shift_success) 
		{
			std::cerr << "\nThere is no merge tree!" << std::endl;
			return false;
		}
	}
	std::cout << std::endl << "FINISHED BOTTOMMOST FILL TREE" << std::endl;
	std::cout<< std::endl << "A bottommost routing Tree is generated!" << std::endl;

	/**********************/
#if defined(LEARNING_ABOUT_TREE)
	unsigned totalNodes = 0;
	for(unsigned int i = 0; i<routingDAG->getNumLayers(); i++)
	{
		Node *layerChecker = routingDAG->findNode(i,0);
		std::cout << "Node " << i << ",0 has in degree: " << layerChecker->getInDegree() << std::endl;			
		Edge *inEdge = layerChecker->getInEdge(0);
		unsigned int j = 0;
		while(layerChecker != NULL)
		{
			layerChecker = routingDAG->findNode(i,j);
			if(layerChecker == NULL)
			{
				continue;
			}
			std::cout << "Node " << i << ", " << j << " thinks it is at ";
			std::cout << layerChecker->getLayer() << ", " << layerChecker->getLayerIndex() << std::endl;
			j++;
		}
		std::cout << "Layer " << i << " has " << j << " total nodes" << std::endl;
		totalNodes +=j;
	}
	std::cout << "The observed total, " << totalNodes << ", should equal " << routingDAG->getNumNodes() << std::endl;
#endif
	/**********************/

	return true;
}

// after the tree is ready, tauten the tree to make it look better
bool UserInterface::tautenTree()
{
	if(!(ctrl_QuadGenerated == QUADRANT_GENERATED && 
		ctrl_WeatherReadIn == WEATHER_READ_IN && 
		ctrl_DemandReadIn == DEMAND_READ_IN &&
		ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED && 
		routingDAG->getStatus()==TREE_GENERATED))
	{
		std::cerr << "Please generate the bottommost tree first before tautening the tree."<<std::endl;
		return false;
	}
	std::cout << "\nWe are tautening the bottommost tree branches now, please wait...";
#if defined(DEBUG_TAUTENING)
	std::cout << "\nBefore tautening" << std::endl;
	for(double rad = .1; rad<25; rad *=2)
	{
		std::cout << "Radius tested: " << rad << '\t';
		routingDAG->areAllNodesFarFromWeather( weatherDataSets, 
												rad, 
												deviationThreshold, 
												nodeEdgeThreshold );
	}
#endif
	bool retval =  routingDAG->generateTautenedTree(weatherDataSets, demandRNPs, deviationThreshold, nodeEdgeThreshold);
#if defined(DEBUG_TAUTENING)
	std::cout << "\nAfter tautening" << std::endl;
	for(double rad = .1; rad<25; rad *=2)
	{
		std::cout << "Radius tested: " << rad << '\t';
		routingDAG->areAllNodesFarFromWeather( weatherDataSets, 
												rad, 
												deviationThreshold, 
												nodeEdgeThreshold );
	}
#endif
	return retval;
}

// After the tree is generated, compute the operational flexibility information for each tree node and edge
void UserInterface::inputOperationalFlexibility()
{
	if(ctrl_RoutingDAGGenerated != ROUTINGDAG_GENERATED || routingDAG->getStatus()==TREE_NOT_GENERATED)
	{
		std::cerr << "\nPlease generate the tree first!"<<std::endl;
		return;
	}
	std::vector<double> radii = inputs.getOperFlex();
	std::sort(radii.begin(), radii.end() );
#if !defined(DO_NOT_CONVERT_LAT_LON_TO_PIXELS)
	for(unsigned int i = 0; i<radii.size(); i++)
	{
		radii[i]/=NMILESPERPIXEL;
	}
#endif
	if( radii[0] < 0 ) 
	{ // valid values, all positive and in increasing order, then move to the next step
		std::cout << "Operational flexbility values are invalid: Verify that all are positive";
		exit(0);
	}
	routingDAG->generateOperFlexPairs(radii, weatherDataSets, deviationThreshold);	// generate the pairs of operational flexibility values	
	ctrl_OperFlexGenerated = OPER_FLEX_GENERATED;
	std::cout << "\nOperational flexibility pairs successfully generated for the tree."<<std::endl;
}

// generate the demand profile for testing
bool UserInterface::inputDemand()
{
	if(ctrl_QuadGenerated==QUADRANT_NOT_GENERATED)
	{
		std::cerr << "\nPlease generate the quadrant first!"<<std::endl;
		return false;
	}
	int numDemands = 0;
	while(true)
	{
		std::cout << "\nPlease input the number of demands (entry nodes):";
		cin>>numDemands;
		if(numDemands>0)
		{
			break;
		}
		else
		{
			std::cout << "\nInvalid Value...";
		}
	}
	std::string tempDemands;
	cin.ignore(1000 ,'\n');													// clear the buffer, ignore the leftover characters
	while(true)
	{
		std::cout << "\nInput the demands one by one in nm (e.g. 2, 2, 1.5, 3, 2.7 means there are 5 entry nodes, and their rnp requirements are 2, 2, 1.5, 3 and 2.7.):";
		getline(cin, tempDemands);
		if(!inputDemandValid(tempDemands, numDemands))		// if the demand format is errorous
		{
			std::cout << "\nInvalid Input...";
		}
		else if (!quadrant->demandFeasible(demandRNPs))		// the new demands are stored in the demandRNPs std::vector but cannot be accomodated by the quadrant
		{
			std::cout << "\nThe quadrant cannot accomodate the new demands...";
		}
		else
		{
			break;
		}
	}
	ctrl_DemandReadIn = DEMAND_READ_IN;
	// when demand or weather data is changed, the routingDAG must be regenerated based on the new demand/weather data
	if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)
	{
		ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;
	}
	return true;
}

// test if the user input demand std::vector is valid input or not, if so, update the demandRNPs std::vector
bool UserInterface::inputDemandValid(std::string &input, int numDemands)
{
	int numCommas = 0;
	int numPoints = 0;
	for(unsigned int i=0; i<input.length(); i++)
	{
		if(input.at(i)==' ')			// eliminate all the spaces
		{
			input.erase(i, 1);
			i--;
		}
	}
	for(unsigned int i=0; i<input.length(); i++)
	{
		char temp = input.at(i);
		if(!((temp>='0' && temp<='9') || temp==' ' || temp==',' || temp=='.'))
		{
			return false;
		}
		if(temp==',')	
		{
			numCommas++;
		}
		if(temp=='.')	
		{
			numPoints++;
		}
		if(i!=input.length()-1 && temp==',' && (input.at(i+1)==',' || input.at(i+1)=='.'))
		{
			return false;							// a comma must be followed by a number
		}
		if(i!=input.length()-1 && temp=='.' && (input.at(i+1)==',' || input.at(i+1)=='.'))
		{
			return false;							// a point must be followed by a number too
		}
	}
	if(numCommas!=numDemands-1)	
	{
		return false;		// there should be numDemands-1 commas in total
	}
	if(numPoints>numDemands)	
	{
		return false;		// there can be at most numDemands floating points
	}
	// after the first lass, the spaces are removed
	char* rnpInputsCStr = new char[input.length()+1];
	for(unsigned int i=0; i<input.length(); i++)
	{
		rnpInputsCStr[i] = input.at(i);
	}
	rnpInputsCStr[input.length()] = '\0';						// almost sure that the input format is correct
	std::vector<double> tempRnps;
	int index = 0;
	int wordStartingIndex = 0;
	int length = strlen(rnpInputsCStr);
	while(index<=length)
	{
		if(rnpInputsCStr[index]!=',' && rnpInputsCStr[index]!='\0')		// parse the std::string by separators
		{
			index++;
		}
		else 
		{
			rnpInputsCStr[index]='\0';
			if(index==wordStartingIndex) 
			{
				return false;			// there is no character between 2 commas, invalid
			}
			int numPoints = 0;
			int floatingPointPosition = 0;
			for(int j=wordStartingIndex; j<index;j++)
			{
				if(rnpInputsCStr[j]=='.')						// for this number, check its doubleing point position
				{
					floatingPointPosition = j;
					numPoints++;
				}
			}
			if(numPoints>1)
			{
				return false;				// one number could have at most 1 doubleing point
			}
			// if there is a floating point, and its position is wrong, return invalid
			if(numPoints==1 && (floatingPointPosition==index-1 || floatingPointPosition==wordStartingIndex))
			{
				return false;
			}
			tempRnps.push_back(atof(&rnpInputsCStr[wordStartingIndex]));
			wordStartingIndex = ++index;
		}
	}
	demandRNPs.clear();											// copy the temporarily stored rnp values into the demand std::vector
	for(unsigned int i=0; i<tempRnps.size(); i++)
	{
		demandRNPs.push_back(tempRnps[i]/NMILESPERPIXEL);
	}
	return true;
}


// read in weather data from a file
bool UserInterface::readWeatherData()
{
	// only when the demand profile is set up and is read in from a file(so that lati/long per pixel are set), are we ready to read in the weather data
	if(ctrl_QuadGenerated==QUADRANT_NOT_GENERATED || ctrl_DemandReadIn == DEMAND_NOT_READ_IN)
	{
		std::cerr << "\nRead In Weather After Setting up the Demand Profile!"<<std::endl;
		return false;
	}
	int totalNumWeatherFiles = inputs.getNumWeatherFiles();
	if(totalNumWeatherFiles<=0) 
	{
		std::cout <<  std::endl << "Invalid weather count..." << std::endl;
		return false;
	}
	std::vector<string> weatherFileDirectories;
	for(int i=0; i<totalNumWeatherFiles; i++)
	{
		std::string tempStr;
		tempStr = inputs.getWeatherFile(i);
		weatherFileDirectories.push_back(tempStr);					// store the directories of each weather file into the std::vector
	}
	/***********************************************************************************************************************************/
	weatherDataSets.clear();
	double rangeMinLati;
	double rangeMinLong;
	double rangeMaxLati;
	double rangeMaxLong;
	// get the range of the weather data by knowing the range of the demand profile, so those unrelated weather
	// data cells are trimed (not read into the memory of the storing std::vector at all)
	demandProfile->getRange(&rangeMinLati, &rangeMinLong, &rangeMaxLati, &rangeMaxLong);
	double minAlt = 10000;
	double maxAlt = 0;			// get the min and std::max altitude of all weather files in order to set the quadrant
	std::cout << "\nReading weather data files now, please wait..." << std::endl;
	/***********************************************************************************************************************************/

	double weatherCellWidth = inputs.getCellWidth();
	for(int i=0; i<totalNumWeatherFiles; i++)
	{
		std::cout << "Parsing weather file: " << weatherFileDirectories[i] << std::endl;
		ifstream is;
		is.open(weatherFileDirectories[i].c_str(), ios::in);	// read in from the weather files
		if( is.is_open() )
		{
			is.close();
			std::string currentFile = weatherFileDirectories[i];
			WeatherData tempWeather;
			if(!tempWeather.readInFileData(currentFile, 
				rangeMinLati-1, 
				rangeMinLong-1, 
				rangeMaxLati+1, 
				rangeMaxLong+1 ) )
			{
				std::cout << "\nWeather not read in successfully..."<<std::endl;
				return false;
			}

			// convert the weather cells to screen OPENGL coordinate system 
			tempWeather.convertLatiLongHeightToXY(centerLati, centerLong, latiPerPixel, longPerPixel, weatherCellWidth);

			weatherDataSets.push_back(tempWeather);		// push the newly read in weather data into the storing std::vector
			minAlt = min(minAlt, (double)tempWeather.getMinAlt());
			maxAlt = std::max(maxAlt, (double)tempWeather.getMaxAlt());
		}
		else
		{
			std::cout << "\nWeather file directory error..."<<std::endl;
			return false;
		}
	} // loop over all weather files
	// weather files are read in, do some tests to make sure that the files are valid
	double totalProbabilityOfWeatherFiles = 0;
	for(unsigned int i=0; i<weatherDataSets.size(); i++)
	{
		totalProbabilityOfWeatherFiles += weatherDataSets[i].getProbability();
		std::cout << i << "   " << weatherDataSets[i].getProbability() << std::endl;
	}
	// the total probability of the weather files is not 1
	if(std::abs(totalProbabilityOfWeatherFiles-1.0) > 0.1)			
	{
		std::cerr <<  "\nThe total probability of the weather data files is not 1."<<std::endl;
		std::cout << "Current Probability: " << totalProbabilityOfWeatherFiles << std::endl;
		weatherDataSets.clear();
		ctrl_WeatherReadIn = WEATHER_NOT_READ_IN;		// the weather is not read in yet
		return false;
	}
	else
	{
		ctrl_WeatherReadIn = WEATHER_READ_IN;
		// when demand or weather data is changed, the routingDAG must be regenerated based on the new demand/weather data
		if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)
		{
			ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;
		}
		quadrant->setiHeight(minAlt);					// if read in successfully, then set the quadrant's altitude information
		quadrant->setoHeight(maxAlt);
	}
	return true;
}

// read in demand profile data from a file
bool UserInterface::readDemandProfile()
{
	/************************************************************************************************/
	if(ctrl_QuadGenerated == QUADRANT_GENERATED)
	{
		std::string demandProfileDirectory;	// the directory of the demand profile
		demandProfileDirectory = inputs.getDemandFile();
		demandProfile->reset();	// a brand new demand instance
		if(demandProfile->readInFile(demandProfileDirectory) )
		{
			ctrl_DemandReadIn = DEMAND_READ_IN;
			// the radius 100 always corresponds to the radius of the transition airspace quadrant, and the conversion from
			// lati/long to screen coordinate system is done by the conversion from 100 to radius in lati/long
			quadrant->setQuadrant(0, 0, inputs.getQuadrantAngle(), inputs.getAngularWidth() , 20, TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL, 0, 0);
			// set up the 4 global variable for drawing purpose (convertion inbetween 2 coordinate systems)
			// the center of the lati/long is always set to be at opengl coordinate (0, 0)
			// and getting the demand profile's starting time and ending time
			demandProfile->getDemandInfo(&centerLati, &centerLong, &latiPerPixel, &longPerPixel, &startTime, &endTime);
			demandRNPs.clear();
			demandProfile->generateDemandVector(demandRNPs, quadrant->getAngle(), quadrant->getAngularWidth(), NMILESPERPIXEL);
			// when demand or weather data is changed, the routingDAG must be regenerated based on the new demand/weather data
			if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)
			{
				ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;
			}
			return true;
		}
		else														// failed to read in successfully
		{
			demandProfile->reset();
		}
	}
	else
	{
		std::cerr << "\nPlease generate a quadrant first!" << std::endl;
	}
	return false;														// the file is NOT read in successfully
}

/**
	\brief after generating the tree, export the tree information into an .xml ascii file
*/
void UserInterface::saveTreeInformation()
{
	if(ctrl_RoutingDAGGenerated == ROUTINGDAG_NOT_GENERATED || routingDAG->getStatus()==TREE_NOT_GENERATED)										
	{
		std::cerr << "\nPlease generate the routing DAG and tree first!"<<std::endl;
		return;
	}
	if(ctrl_OperFlexGenerated == OPER_FLEX_NOT_GENERATED)
	{
		std::cerr << "\nPlease generate Operational Flexity Pairs First!"<<std::endl;
		return;
	}
	routingDAG->outputTreeInformation(centerLati, centerLong, latiPerPixel, longPerPixel, startTime, endTime, outputFileName, weatherDataSets, deviationThreshold );
	std::cout << "latiPerPixel: " << latiPerPixel << "; centerLati: " << centerLati << std::endl;
	std::cout << "longPerPixel: " << longPerPixel << "; centerLong: " << centerLong << std::endl;
	std::cout<< "\nTree Information successfully written to file.\n";
}

