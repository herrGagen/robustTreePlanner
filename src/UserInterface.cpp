#include "UserInterface.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
/*****************************************/
// header files from within the project
#include "NodeAndEdge.h"

using namespace std;	

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

// used to generate a brand new tree generating instance, everything is restalled to original value
void UserInterface::reset()
{
  // deleting the existing objects
  if(ctrl_QuadGenerated == QUADRANT_GENERATED)			delete quadrant;
  if(ctrl_WeatherReadIn == WEATHER_READ_IN)				delete demandProfile;
  if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)	delete routingDAG;
  weatherData.clear();
  demandRNPs.clear();
  resetHelper();											// set all the ctrl variables to default and generate new objects
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
  deviationThreshold = (float) 0.7;
  nodeEdgeThreshold = (float) 0.8;								// default values of the thresholds
}


/****************************************************************************************/
// RTPrototype Functions

// The project starts excuting here, called by the main function
void UserInterface::ProgramBegins(std::string inputFile)
{
  ifstream in_stream;
  string line;

  std::cout << "Input file: " << inputFile.c_str() << endl;
  in_stream.open(inputFile.c_str());
  currentInput = 0;
  while (!in_stream.eof() ) {
    getline(in_stream, line);
    for(int j = line.length()-1; j>=0; j--)
    {
      if(line[j] == '\n' || line[j] == '\r')
      {
        line[j] = '\0';
      }
      else
      {
        break;
      }
    }
    allInputs[currentInput] = line;
    ++currentInput;
  }
  in_stream.close();
  currentInput = 0;

  char userInput = '0';	     // the user input will always be a number, deciding what to do next
  bool startOver = false;    // when the user choose to do everything again
  std::cout<<"Robust Tree Generator:\nThe software is used to generate robust trees given the demand profile and weather data information."<<endl;

  if(startOver) // when doing everything again, reset everything first
  {
    reset();
    startOver = false;
  }
  if(!readDemandProfile()) // at the start, read in demand profile
  {
    std::cout<<"\nFailed to read in the demand profile, double check the demand path listed."<<endl;
    return;
  }
  std::cout<<"Demand profile successfully read in." << endl;
  if(!readWeatherData())
  {
    std::cout<<"\nFailed to successfully read in the weather data, double check the weather directories."<<endl;
    return;
  }
  std::cout<<"\nWeather files are succesfully read in!"<<endl;
  /*  (*quadrant).setQuadrant( quadrant->getcX(),
  quadrant->getcY(),
  ::atof(allInputs[currentInput++].c_str()),
  quadrant->getiRadius(),
  quadrant->getoRadius(),
  quadrant->getiHeight(),
  quadrant->getoHeight() ); */
  (*quadrant).setAngle( ::atof(allInputs[currentInput++].c_str() ) );
  //   std::cout << "Current LINE FOR allInputs: " << currentInput << "  " << allInputs[currentInput] << endl;
  std::cout << "Current angle: " << quadrant->getAngle() << endl;
  std::cout << "Generating tree." << endl;
  generateTree();
  std::cout << "Tautening tree." << endl;
  tautenTree();
  std::cout << "Doing operational Flexibility stuff." << endl;
  inputOperationalFlexibility();
  std::cout << "Saving tree information." << endl;
  saveTreeInformation();

}

// print the information of the quadrant and the demand rnps for each entry node
void UserInterface::printQuadrantAndDemandInfo()
{
  if(ctrl_QuadGenerated == QUADRANT_GENERATED)  
  {
    // first print the current quadrant information and demand information
    std::cout<<"\nThe center of the quadrant is currently at lati/long ("<<setprecision(5)<<centerLati+quadrant->getcX()*latiPerPixel<<", ";
    std::cout<<setprecision(5)<<centerLong+quadrant->getcY()*longPerPixel<<").";
    std::cout<<"\nThe inner radius (in nm) is "<<setprecision(5)<<quadrant->getiRadius()*NMILESPERPIXEL<<", the outer radius is ";
    std::cout<<setprecision(5)<<quadrant->getoRadius()*NMILESPERPIXEL<<".";
    std::cout<<"\nThe inner altitude (in feet) is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getiHeight()*ALTITUDE_PER_PIXEL<<", ";
    std::cout<<"the outer altitude is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getoHeight()*ALTITUDE_PER_PIXEL<<".";
    std::cout<<"\nThe angle of the right boundary of the quadrant is (relative to +x axis) "<<setprecision(5)<<quadrant->getAngle()*180/PI<<" degrees."<<endl;
    // then print the demand information
    std::cout<<"\nCurrently there are "<<demandRNPs.size()<<" entry points evenly spaced on the outer boundary of the quadrant.";
    if(!demandRNPs.empty())
    {
      std::cout<<"\nthe rnp requirements are (in nm, from right to left) ";
      for(unsigned i=0; i<demandRNPs.size(); i++)
        std::cout<<setprecision(5)<<demandRNPs[i]*NMILESPERPIXEL<<"  ";
    }
  }
  else
    std::cout<<"\nThe quadrant is not generated yet!"<<endl;
}

// display the quadrant infomation and prompt the users if they want to edit the information
bool UserInterface::editQuadrant()
{
  if(ctrl_QuadGenerated == QUADRANT_GENERATED)  
  {
    std::cout<<"Please input the new quadrant information";
    double cx, cy, angle, ir, ora, ih, oh;
    std::cout<<"\nThe latitude of the center of the quadrant is (enter 0 to skip):";
    cin>>cx;
    if(cx!=0)	cx = (cx-centerLati) / latiPerPixel;
    std::cout<<"\nThe longitude of the center of the quadrant is (enter 0 to skip):";
    cin>>cy;
    if(cy!=0)	cy = (cy-centerLong) / longPerPixel;
    std::cout<<"\nThe angle of the right side of the quadrant is (0 to 360 degrees, relative to the positive x axis, enter 0 to skip):";
    cin>> angle;
    if(angle!=0) angle = angle*PI/180;
    std::cout<<"\nThe inner radius of the quadrant is (in nm, enter 0 to skip):";
    cin>>ir;
    if(ir!=0)	ir = abs(ir)/NMILESPERPIXEL;
    std::cout<<"\nThe outer radius of the quadrant is (in nm, enter 0 to skip):";
    cin>>ora;
    if(ora!=0)	ora = abs(ora)/NMILESPERPIXEL;
    std::cout<<"\nThe inner altitude of the quadrant is (in feet, enter 0 to skip):";
    cin>>ih;
    if(ih!=0)	ih = (ih-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
    std::cout<<"\nThe outer altitude of the quadrant is (in feet, enter 0 to skip):";
    cin>>oh;
    if(oh!=0)	oh = (oh-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
    // set the new information gotten from the dialog to the quadrant object
    quadrant->setQuadrant(cx, cy, angle, ir, ora, ih, oh);
    printQuadrantAndDemandInfo();										// after editing, print the new quadrant information out
    if(ctrl_DemandReadIn == DEMAND_READ_IN && !quadrant->demandFeasible(demandRNPs))
    {
      demandRNPs.clear();
      ctrl_DemandReadIn = DEMAND_NOT_READ_IN;	// if the new quadrant cannot accomodate the demands, then set demand NOT read in
      std::cout<<"\nThe new quadrant cannot accommondate the rnp requirements of the entry nodes, the demand info are therefore restored."<<endl;
      return false;
    }
    // quadrant is changed, then if a tree or routing DAG was generated, they are considered outdated
    ctrl_RoutingDAGGenerated = !ROUTINGDAG_GENERATED;	
  }
  else
    std::cout<<"\nThe quadrant is not generated yet!"<<endl;
  return true;
}

/*****************************************************************************************************************************************/

/** 
\brief Generates a tree
*/
bool UserInterface::generateTree()
{
  // only when the quadrant is generated and weather data is read in can we generate a tree
  if(ctrl_QuadGenerated == QUADRANT_GENERATED && ctrl_WeatherReadIn == WEATHER_READ_IN && ctrl_DemandReadIn == DEMAND_READ_IN)
  {
    double minDistBetweenMergeNodes = 5/NMILESPERPIXEL; // WILLXYZ make this read from config file

    std::cout << "Reading in deviation threshold." << endl;
    deviationThreshold = ::atof(allInputs[currentInput++].c_str());

    std::cout << "Reading in nodeEdgeThreshold." << endl;
    nodeEdgeThreshold = ::atof(allInputs[currentInput++].c_str());
    routingDAG->setminimumDistanceBetweenMergingNodes(minDistBetweenMergeNodes);
    std::cout<< endl << "Generating a bottommost merge tree, please wait..." << endl;
    /************************************************************************************************/
    ctrl_OperFlexGenerated = OPER_FLEX_NOT_GENERATED;				// when generating a new tree, the Oper-Flex pairs need to be generated again
    // first, generate the entry and fix nodes, then the internal nodes
    routingDAG->reset();											// a brand new routing instance
    std::cout << endl << "Finished resetting edges." << endl;

    std::cout << "Current Input value: " << allInputs[currentInput] << endl;
    double quadrantAngleOffset = ::atof(allInputs[currentInput++].c_str());
    double lane_width = ::atof(allInputs[currentInput++].c_str());
    if (lane_width > 0) {
      for (int i = 0; i < demandRNPs.size(); i++) {
        if (demandRNPs[i] > 0) {
          demandRNPs[i] = lane_width;
        }
      }
    }
    std::cout << "Demand RNPS: ";
    for(int i = 0; i < demandRNPs.size(); i++) {
      std::cout << demandRNPs[i] << " ";
    }
    std::cout << endl;

    std::cout << "Max Fix Nodes: " << allInputs[currentInput] << endl;
    int maxFixNodes = ::atof(allInputs[currentInput++].c_str());
    std::cout << "START GENERATING DAG" << endl;

    int demand_shift = ::atoi(allInputs[currentInput++].c_str());
    int demand_drop  = ::atoi(allInputs[currentInput++].c_str());

    if(quadrant->generateDAG(demandRNPs, demandRNPs.size(), deviationThreshold, nodeEdgeThreshold, weatherData, routingDAG, quadrantAngleOffset, maxFixNodes))
    {
      std::cout << "FINISHED DAG" << endl;
      std::cout << endl << "Generating edge set..." << endl;
      routingDAG->generateEdgeSet();								// generate the edges in the searching DAG
      std::cout << endl << "Edges (Routing DAG) generated." << endl;
      ctrl_RoutingDAGGenerated = ROUTINGDAG_GENERATED;
      std::cout << endl << "Generating Tree..." << endl;
      // generate the tree here
      std::cout << "START GENERATING TREE" << endl << "Demand RNPs: ";
      for (int i = 0; i < demandRNPs.size(); i++) {
        std::cout << demandRNPs[i] << "  ";
      }
      std::cout << endl;

      if(!routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold))
      {
        // Will added some code 04/2013 to add a demand shifting scheme to hopefully reduce tree generation failure

        int demand_shift_index = 0;
        bool demand_shift_success = false;
        if (demand_shift == 1) {
          std::cout << endl << "Beginning demand shifting." << endl;
          std::cout << "Variables: demand_shift " << demand_shift << "; demand drop " << demand_drop << "; demand_shift_index " << demand_shift_index << "; demand_shift_success " << demand_shift_success << endl;
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
          while (demand_shift_index < demandRNPs.size() && !demand_shift_success ) {

            if (demand_shift_index > 0 && demandRNPs[demand_shift_index - 1] == 0 && !demand_shift_success) {
              // Swap the values of demandRNPs at the current index and the one below
              demandRNPs[demand_shift_index - 1] = demandRNPs[demand_shift_index];
              demandRNPs[demand_shift_index] = 0;

              // If the tree is generated successfully, exit this block gracefully
              if (routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold)) {
                demand_shift_success = true;
                std::cout << endl << "A downward swap generated a tree on the following index: " << demand_shift_index << endl;
              }

              // Swap them back
              demandRNPs[demand_shift_index] = demandRNPs[demand_shift_index - 1];
              demandRNPs[demand_shift_index - 1] = 0;

            }
            if (demand_shift_index + 1 < demandRNPs.size() && demandRNPs[demand_shift_index + 1] == 0 && !demand_shift_success) {
              // Swap the values of demandRNPs at the current index and the one above
              demandRNPs[demand_shift_index + 1] = demandRNPs[demand_shift_index];
              demandRNPs[demand_shift_index] = 0;

              // If the tree is generated successfully, exit this block gracefully
              if (routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold)) {
                demand_shift_success = true;
                std::cout << endl << "An upward swap generated a tree on the following index: " << demand_shift_index << endl;
              }

              // Swap them back
              demandRNPs[demand_shift_index] = demandRNPs[demand_shift_index + 1];
              demandRNPs[demand_shift_index + 1] = 0;
            }
            demand_shift_index++;
          }
        }
        if (demand_drop >= 1 && !demand_shift_success) {
          /* Begin phase (2) of the demand shifting, aka demand_drop */
          std::cout << endl <<  "Demand shifting unsuccessful, beginning demand dropping." << endl;

          int comb_i = 0; // an index of combinations
          int k = 0; // k is the last positive element
          float temp_demands[4]; // This 4 needs to be changed if combinations is ever changed -- HARDCODED
          // the below 255 is the length of ``combinations`` -- HARDCODED
          while (k < demand_drop && comb_i < 255 && !demand_shift_success) {

            for (int j = 0; j < 4; ++j) {
              temp_demands[j] = 0;
            }

            // std::cout << endl;
            for (int j = 0; j < 4; ++j) {
              if (demandRNPs[combinations[comb_i][j]] > 0) {
                temp_demands[j] = demandRNPs[combinations[comb_i][j]];
                demandRNPs[combinations[comb_i][j]] = 0;
              }
              else {
                // We could break out of the loop, but then we also need to reset demandRNPs values
                // break;
              }
              // HARDCODED VALUE
              if (j == 3) {
                if (routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold)) {
                  demand_shift_success = true;
                  std::cout << endl << "A demand drop generated a tree on the following indices: ";
                  for(int i = 0; i <= k; ++i) {
                    std::cout << combinations[comb_i][i] << " ";
                  }
                  std::cout << endl;

                  std::cout << endl << "Demand RNPs after the drop:" << endl;
                  for(int i = 0; i < demandRNPs.size(); i++) {
                    std::cout << demandRNPs[i] << " ";
                  }
                  std::cout << endl;
                }
                else {
                  // reset all of the demandRNPs values -- HARDCODED
                  /*
                  std::cout << "Demand RNPS: ";
                  for(int i = 0; i < demandRNPs.size(); i++) {
                  std::cout << demandRNPs[i] << " ";
                  }
                  */
                  for (int j = 0; j < 4; ++j) {
                    demandRNPs[combinations[comb_i][j]] = temp_demands[j];
                    if (demandRNPs[combinations[comb_i][j]] < 0.00000001) {
                      demandRNPs[combinations[comb_i][j]] = 0;
                    }
                  }
                }
              }
            }
            // Find the new k value, which is the last positive number in the array -- HARDCODED
            for (int j = 0; j < 4; ++j) {
              if (combinations[comb_i][j] < 0) {
                k = j - 1;
                break;

                if (j == 3) { k = 3; } // j would go on to being 4, but can't because of the loop condition
              }
            }
            ++comb_i;
          }
          /*
          float temp_demand;
          while (demand_shift_index < demandRNPs.size() && !demand_shift_success) {
          if (demandRNPs[demand_shift_index] > 0) {
          temp_demand = demandRNPs[demand_shift_index];
          demandRNPs[demand_shift_index] = 0;

          if (routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold)) {
          demand_shift_success = true;
          std::cout << "A demand drop generated a tree on the following index: " << demand_shift_index << endl;
          }

          demandRNPs[demand_shift_index] = temp_demand;
          }
          demand_shift_index++;
          }
          */
        }

        if (!demand_shift_success) {
          cerr<< endl << "There Does NOT Exist A Merge Tree!"<<endl;
          exit(0);
          return false;
        }
      }
      std::cout << endl << "FINISHED BOTTOMMOST FILL TREE" << endl;
      std::cout<< endl << "A bottommost routing Tree is generated!" << endl;
      return true;
    }	// an error message will pop up if failed to generate the DAG
    else { std::cout << "Failed to generate the DAG."; }
  }
  // else, then the weather data and demand profile have to be read in first
  else cerr<<"\nPlease read in or generate the demand profile and weather data first."<<endl;
  return false;
}

// after the tree is ready, tauten the tree to make it look better
bool UserInterface::tautenTree()
{
  if(!(ctrl_QuadGenerated == QUADRANT_GENERATED && ctrl_WeatherReadIn == WEATHER_READ_IN && ctrl_DemandReadIn == DEMAND_READ_IN &&
    ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED && routingDAG->getStatus()==TREE_GENERATED))
  {
    cerr<<"Please generate the bottommost tree first before tautening the tree."<<endl;
    return false;
  }
  std::cout<<"\nWe are tautening the bottommost tree branches now, please wait...";
  return routingDAG->generateTautenedTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold);
}

// After the tree is generated, compute the operational flexibility information for each tree node and edge
void UserInterface::inputOperationalFlexibility()
{
  if(ctrl_RoutingDAGGenerated != ROUTINGDAG_GENERATED || routingDAG->getStatus()==TREE_NOT_GENERATED)
  {
    cerr<<"\nPlease generate the tree first!"<<endl;
    return;
  }

  // XYZ WILL changed this on 2013-04-30
  // Cleaning code up and making the number of operational flexibility parameters variable
  int numOperFlexParams = ::atoi(allInputs[currentInput++].c_str());

  float* radii = new float[numOperFlexParams];
  for (int i = 0; i < numOperFlexParams; ++i) {
    radii[i] = ::atof(allInputs[currentInput++].c_str());
  }
  routingDAG->generateOperFlexPairs(radii, numOperFlexParams, weatherData, deviationThreshold);	// generate the pairs of operational flexibility values	
  delete []radii;
  ctrl_OperFlexGenerated = OPER_FLEX_GENERATED;
  std::cout<<"\nOperational flexibility pairs successfully generated for the tree."<<endl;
}

// generate the demand profile for testing
bool UserInterface::inputDemand()
{
  if(ctrl_QuadGenerated==QUADRANT_NOT_GENERATED)
  {
    cerr<<"\nPlease generate the quadrant first!"<<endl;
    return false;
  }
  int numDemands = 0;
  while(true)
  {
    std::cout<<"\nPlease input the number of demands (entry nodes):";
    cin>>numDemands;
    if(numDemands>0)
      break;
    else
      std::cout<<"\nInvalid Value...";
  }
  string tempDemands;
  cin.ignore(1000 ,'\n');													// clear the buffer, ignore the leftover characters
  while(true)
  {
    std::cout<<"\nInput the demands one by one in nm (e.g. 2, 2, 1.5, 3, 2.7 means there are 5 entry nodes, and their rnp requirements are 2, 2, 1.5, 3 and 2.7.):";
    getline(cin, tempDemands);
    if(!inputDemandValid(tempDemands, numDemands))		// if the demand format is errorous
      std::cout<<"\nInvalid Input...";
    else if (!quadrant->demandFeasible(demandRNPs))		// the new demands are stored in the demandRNPs vector but cannot be accomodated by the quadrant
      std::cout<<"\nThe quadrant cannot accomodate the new demands...";
    else
      break;
  }
  ctrl_DemandReadIn = DEMAND_READ_IN;
  // when demand or weather data is changed, the routingDAG must be regenerated based on the new demand/weather data
  if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)
    ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;
  return true;
}

// test if the user input demand vector is valid input or not, if so, update the demandRNPs vector
bool UserInterface::inputDemandValid(string &input, int numDemands)
{
  int numCommas = 0;
  int numPoints = 0;
  for(int i=0; i<input.length(); i++)
  {
    if(input.at(i)==' ')			// eliminate all the spaces
    {
      input.erase(i, 1);
      i--;
    }
  }
  for(int i=0; i<input.length(); i++)
  {
    char temp = input.at(i);
    if(!((temp>='0' && temp<='9') || temp==' ' || temp==',' || temp=='.'))
      return false;
    if(temp==',')	numCommas++;
    if(temp=='.')	numPoints++;
    if(i!=input.length()-1 && temp==',' && (input.at(i+1)==',' || input.at(i+1)=='.'))
      return false;							// a comma must be followed by a number
    if(i!=input.length()-1 && temp=='.' && (input.at(i+1)==',' || input.at(i+1)=='.'))
      return false;							// a point must be followed by a number too
  }
  if(numCommas!=numDemands-1)	return false;		// there should be numDemands-1 commas in total
  if(numPoints>numDemands)	return false;		// there can be at most numDemands floating points
  // after the first lass, the spaces are removed
  char* rnpInputsCStr = new char[input.length()+1];
  for(int i=0; i<input.length(); i++)
    rnpInputsCStr[i] = input.at(i);
  rnpInputsCStr[input.length()] = '\0';						// almost sure that the input format is correct
  vector<float> tempRnps;
  int index = 0;
  int wordStartingIndex = 0;
  int length = strlen(rnpInputsCStr);
  while(index<=length)
  {
    if(rnpInputsCStr[index]!=',' && rnpInputsCStr[index]!='\0')		// parse the string by separators
      index++;
    else 
    {
      rnpInputsCStr[index]='\0';
      if(index==wordStartingIndex) return false;			// there is no character between 2 commas, invalid
      int numPoints = 0;
      int floatingPointPosition = 0;
      for(int j=wordStartingIndex; j<index;j++)
      {
        if(rnpInputsCStr[j]=='.')						// for this number, check its floating point position
        {
          floatingPointPosition = j;
          numPoints++;
        }
      }
      if(numPoints>1)			return false;				// one number could have at most 1 floating point
      // if there is a floating point, and its position is wrong, return invalid
      if(numPoints==1 && (floatingPointPosition==index-1 || floatingPointPosition==wordStartingIndex))
        return false;
      tempRnps.push_back(atof(&rnpInputsCStr[wordStartingIndex]));
      wordStartingIndex = ++index;
    }
  }
  demandRNPs.clear();											// copy the temporarily stored rnp values into the demand vector
  for(int i=0; i<tempRnps.size(); i++)
    demandRNPs.push_back(tempRnps[i]/NMILESPERPIXEL);
  return true;
}


// read in weather data from a file
bool UserInterface::readWeatherData()
{
  // only when the demand profile is set up and is read in from a file(so that lati/long per pixel are set), are we ready to read in the weather data
  if(ctrl_QuadGenerated==QUADRANT_NOT_GENERATED || ctrl_DemandReadIn == DEMAND_NOT_READ_IN)
  {
    cerr<<"\nRead In Weather After Setting up the Demand Profile!"<<endl;
    return false;
  }
  int totalNumWeatherFiles;
  istringstream is(allInputs[currentInput++]);
  is >> totalNumWeatherFiles;
  if(totalNumWeatherFiles<=0) {
    std::cout<< endl << "Invalid weather count..." << endl;
    return false;
  }
  vector<string> weatherFileDirectories;
  for(int i=0; i<totalNumWeatherFiles; i++)
  {
    string tempStr;
    tempStr = allInputs[currentInput++];
    weatherFileDirectories.push_back(tempStr);					// store the directories of each weather file into the vector
  }
  /***********************************************************************************************************************************/
  weatherData.clear();
  double rangeMinLati, rangeMinLong, rangeMaxLati, rangeMaxLong;
  // get the range of the weather data by knowing the range of the demand profile, so those unrelated weather
  // data cells are trimed (not read into the memory of the storing vector at all)
  demandProfile->getRange(&rangeMinLati, &rangeMinLong, &rangeMaxLati, &rangeMaxLong);
  float minAlt, maxAlt;	minAlt = 10000;		maxAlt = 0;			// get the min and max altitude of all weather files in order to set the quadrant
  std::cout<<"\nReading weather data files now, please wait..."<<endl;
  /***********************************************************************************************************************************/

  double weatherCellWidth = ::atof(allInputs[currentInput++].c_str());
  for(int i=0; i<totalNumWeatherFiles; i++)
  {
    std::cout << "Parsing weather file: " << weatherFileDirectories[i] << std::endl;
    ifstream is(weatherFileDirectories[i].c_str(), ios::in);	// read in from the weather files
    if(true) //is.is_open())
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
        std::cout<<"\nWeather not read in successfully..."<<endl;
        return false;
      }

      // You had something about onverting lat long to pixel values, but used uninitialized values.
      tempWeather.convertLatiLongHeightToXY(centerLati, centerLong, latiPerPixel, longPerPixel, weatherCellWidth);

      weatherData.push_back(tempWeather);		// push the newly read in weather data into the storing vector
      minAlt = min(minAlt, (float)tempWeather.getMinAlt());
      maxAlt = max(maxAlt, (float)tempWeather.getMaxAlt());
    }
    else
    {
      std::cout<<"\nWeather file directory error..."<<endl;
      return false;
    }
  } // loop over all weather files
  // weather files are read in, do some tests to make sure that the files are valid
  float totalProbabilityOfWeatherFiles = 0;
  for(int i=0; i<weatherData.size(); i++)
  {
    totalProbabilityOfWeatherFiles += weatherData[i].getProbability();
    std::cout << i << "   " << weatherData[i].getProbability() << endl;
  }
  // the total probability of the weather files is not 1
  if(abs(totalProbabilityOfWeatherFiles-1.0)>0.1)			
  {
    cerr<< "\nThe total probability of the weather data files is not 1."<<endl;
    std::cout << "Current Probability: " << totalProbabilityOfWeatherFiles << endl;
    weatherData.clear();
    ctrl_WeatherReadIn = WEATHER_NOT_READ_IN;		// the weather is not read in yet
    return false;
  }
  else
  {
    ctrl_WeatherReadIn = WEATHER_READ_IN;
    // when demand or weather data is changed, the routingDAG must be regenerated based on the new demand/weather data
    if(ctrl_RoutingDAGGenerated == ROUTINGDAG_GENERATED)
      ctrl_RoutingDAGGenerated = ROUTINGDAG_NOT_GENERATED;
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
    string demandProfileDirectory;	// the directory of the demand profile
    demandProfileDirectory = allInputs[currentInput++];
    demandProfile->reset();	// a brand new demand instance
    if(demandProfile->readInFile(demandProfileDirectory) )
    {
      ctrl_DemandReadIn = DEMAND_READ_IN;
      // the radius 100 always corresponds to the radius of the transition airspace quadrant, and the conversion from
      // lati/long to screen coordinate system is done by the conversion from 100 to radius in lati/long
      quadrant->setQuadrant(0, 0, 10*PI/180, 20, TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL, 0, 0);
      // set up the 4 global variable for drawing purpose (convertion inbetween 2 coordinate systems)
      // the center of the lati/long is always set to be at opengl coordinate (0, 0)
      // and getting the demand profile's starting time and ending time
      demandProfile->getDemandInfo(&centerLati, &centerLong, &latiPerPixel, &longPerPixel, &startTime, &endTime);
      demandRNPs.clear();
      demandProfile->generateDemandVector(demandRNPs, quadrant->getAngle(), quadrant->getAngle()+PI/2, NMILESPERPIXEL);
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
    cerr<<"\nPlease generate a quadrant first!"<<endl;
  //  std::cout << ctrl_QuadGenerated << endl << QUADRANT_GENERATED << endl;
  return false;														// the file is NOT read in successfully
}

// after generating the tree, export the tree information into an .xml ascii file
void UserInterface::saveTreeInformation()
{
  if(ctrl_RoutingDAGGenerated == ROUTINGDAG_NOT_GENERATED || routingDAG->getStatus()==TREE_NOT_GENERATED)										
  {
    cerr<<"\nPlease generate the routing DAG and tree first!"<<endl;
    return;
  }
  if(ctrl_OperFlexGenerated == OPER_FLEX_NOT_GENERATED)
  {
    cerr<<"\nPlease generate Operational Flexity Pairs First!"<<endl;
    return;
  }
  routingDAG->outputTreeInformation(centerLati, centerLong, latiPerPixel, longPerPixel, startTime, endTime, allInputs[currentInput++]);
  std::cout << "latiPerPixel: " << latiPerPixel << "; centerLati: " << centerLati << endl;
  std::cout << "longPerPixel: " << longPerPixel << "; centerLong: " << centerLong << endl;
  std::cout<< "\nTree Information successfully written to file.\n";
}

