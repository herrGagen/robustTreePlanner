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
  cout<<"Robust Tree Generator:\nThe software is used to generate robust trees given the demand profile and weather data information."<<endl;

  if(startOver) // when doing everything again, reset everything first
  {
    reset();
    startOver = false;
  }
  if(!readDemandProfile()) // at the start, read in demand profile
  {
    cout<<"\nFailed to read in the demand profile, double check the demand path listed."<<endl;
    return;
  }
  cout<<"Demand profile successfully read in." << endl;
  if(!readWeatherData())
  {
    cout<<"\nFailed to successfully read in the weather data, double check the weather directories."<<endl;
    return;
  }
  cout<<"\nWeather files are succesfully read in!"<<endl;
  /*  (*quadrant).setQuadrant( quadrant->getcX(),
  quadrant->getcY(),
  ::atof(allInputs[currentInput++].c_str()),
  quadrant->getiRadius(),
  quadrant->getoRadius(),
  quadrant->getiHeight(),
  quadrant->getoHeight() ); */
  (*quadrant).setAngle( ::atof(allInputs[currentInput++].c_str() ) );
  //   cout << "Current LINE FOR allInputs: " << currentInput << "  " << allInputs[currentInput] << endl;
  cout << "Current angle: " << quadrant->getAngle() << endl;
  cout << "Generating tree." << endl;
  generateTree();
  cout << "Tautening tree." << endl;
  tautenTree();
  cout << "Doing operational Flexibility stuff." << endl;
  inputOperationalFlexibility();
  cout << "Saving tree information." << endl;
  saveTreeInformation();

}

// print the information of the quadrant and the demand rnps for each entry node
void UserInterface::printQuadrantAndDemandInfo()
{
  if(ctrl_QuadGenerated == QUADRANT_GENERATED)  
  {
    // first print the current quadrant information and demand information
    cout<<"\nThe center of the quadrant is currently at lati/long ("<<setprecision(5)<<centerLati+quadrant->getcX()*latiPerPixel<<", ";
    cout<<setprecision(5)<<centerLong+quadrant->getcY()*longPerPixel<<").";
    cout<<"\nThe inner radius (in nm) is "<<setprecision(5)<<quadrant->getiRadius()*NMILESPERPIXEL<<", the outer radius is ";
    cout<<setprecision(5)<<quadrant->getoRadius()*NMILESPERPIXEL<<".";
    cout<<"\nThe inner altitude (in feet) is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getiHeight()*ALTITUDE_PER_PIXEL<<", ";
    cout<<"the outer altitude is "<<setprecision(5)<<ALTITUDE_AT_BASE_PLANE+quadrant->getoHeight()*ALTITUDE_PER_PIXEL<<".";
    cout<<"\nThe angle of the right boundary of the quadrant is (relative to +x axis) "<<setprecision(5)<<quadrant->getAngle()*180/PI<<" degrees."<<endl;
    // then print the demand information
    cout<<"\nCurrently there are "<<demandRNPs.size()<<" entry points evenly spaced on the outer boundary of the quadrant.";
    if(!demandRNPs.empty())
    {
      cout<<"\nthe rnp requirements are (in nm, from right to left) ";
      for(unsigned i=0; i<demandRNPs.size(); i++)
        cout<<setprecision(5)<<demandRNPs[i]*NMILESPERPIXEL<<"  ";
    }
  }
  else
    cout<<"\nThe quadrant is not generated yet!"<<endl;
}

// display the quadrant infomation and prompt the users if they want to edit the information
bool UserInterface::editQuadrant()
{
  if(ctrl_QuadGenerated == QUADRANT_GENERATED)  
  {
    cout<<"Please input the new quadrant information";
    double cx, cy, angle, ir, ora, ih, oh;
    cout<<"\nThe latitude of the center of the quadrant is (enter 0 to skip):";
    cin>>cx;
    if(cx!=0)	cx = (cx-centerLati) / latiPerPixel;
    cout<<"\nThe longitude of the center of the quadrant is (enter 0 to skip):";
    cin>>cy;
    if(cy!=0)	cy = (cy-centerLong) / longPerPixel;
    cout<<"\nThe angle of the right side of the quadrant is (0 to 360 degrees, relative to the positive x axis, enter 0 to skip):";
    cin>> angle;
    if(angle!=0) angle = angle*PI/180;
    cout<<"\nThe inner radius of the quadrant is (in nm, enter 0 to skip):";
    cin>>ir;
    if(ir!=0)	ir = abs(ir)/NMILESPERPIXEL;
    cout<<"\nThe outer radius of the quadrant is (in nm, enter 0 to skip):";
    cin>>ora;
    if(ora!=0)	ora = abs(ora)/NMILESPERPIXEL;
    cout<<"\nThe inner altitude of the quadrant is (in feet, enter 0 to skip):";
    cin>>ih;
    if(ih!=0)	ih = (ih-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
    cout<<"\nThe outer altitude of the quadrant is (in feet, enter 0 to skip):";
    cin>>oh;
    if(oh!=0)	oh = (oh-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
    // set the new information gotten from the dialog to the quadrant object
    quadrant->setQuadrant(cx, cy, angle, ir, ora, ih, oh);
    printQuadrantAndDemandInfo();										// after editing, print the new quadrant information out
    if(ctrl_DemandReadIn == DEMAND_READ_IN && !quadrant->demandFeasible(demandRNPs))
    {
      demandRNPs.clear();
      ctrl_DemandReadIn = DEMAND_NOT_READ_IN;	// if the new quadrant cannot accomodate the demands, then set demand NOT read in
      cout<<"\nThe new quadrant cannot accommondate the rnp requirements of the entry nodes, the demand info are therefore restored."<<endl;
      return false;
    }
    // quadrant is changed, then if a tree or routing DAG was generated, they are considered outdated
    ctrl_RoutingDAGGenerated = !ROUTINGDAG_GENERATED;	
  }
  else
    cout<<"\nThe quadrant is not generated yet!"<<endl;
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

    cout << "Reading in deviation threshold." << endl;
    deviationThreshold = ::atof(allInputs[currentInput++].c_str());

    cout << "Reading in nodeEdgeThreshold." << endl;
    nodeEdgeThreshold = ::atof(allInputs[currentInput++].c_str());
    routingDAG->setminimumDistanceBetweenMergingNodes(minDistBetweenMergeNodes);
    cout<< endl << "Generating a bottommost merge tree, please wait..." << endl;
    /************************************************************************************************/
    ctrl_OperFlexGenerated = OPER_FLEX_NOT_GENERATED;				// when generating a new tree, the Oper-Flex pairs need to be generated again
    // first, generate the entry and fix nodes, then the internal nodes
    routingDAG->reset();											// a brand new routing instance
    cout << endl << "Finished resetting edges." << endl;

    cout << "Current Input value: " << allInputs[currentInput] << endl;
    double quadrantAngleOffset = ::atof(allInputs[currentInput++].c_str());
    double lane_width = ::atof(allInputs[currentInput++].c_str());
    if (lane_width > 0) {
      for (int i = 0; i < demandRNPs.size(); i++) {
        if (demandRNPs[i] > 0) {
          demandRNPs[i] = lane_width;
        }
      }
    }
    cout << "Demand RNPS: ";
    for(int i = 0; i < demandRNPs.size(); i++) {
      cout << demandRNPs[i] << " ";
    }
    cout << endl;

    cout << "START GENERATING DAG" << endl;


    if(quadrant->generateDAG(demandRNPs, demandRNPs.size(), deviationThreshold, nodeEdgeThreshold, weatherData, routingDAG, quadrantAngleOffset))
    {
      cout << "FINISHED DAG" << endl;
      cout << endl << "Generating edge set..." << endl;
      routingDAG->generateEdgeSet();								// generate the edges in the searching DAG
      cout << endl << "Edges (Routing DAG) generated." << endl;
      ctrl_RoutingDAGGenerated = ROUTINGDAG_GENERATED;
      cout << endl << "Generating Tree..." << endl;
      // generate the tree here
      cout << "START GENERATING TREE" << endl;
      if(!routingDAG->generateTree(weatherData, demandRNPs, deviationThreshold, nodeEdgeThreshold))
      {
        cerr<< endl << "There Does NOT Exist A Merge Tree!"<<endl;
        return false;
      }
      else 
      {
        cout << endl << "FINISHED BOTTOMMOST FILL TREE" << endl;
        cout<< endl << "A bottommost routing Tree is generated!" << endl;
        return true;
      }
    }	// an error message will pop up if failed to generate the DAG
    else { cout << "Failed to generate the DAG."; }
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
  cout<<"\nWe are tautening the bottommost tree branches now, please wait...";
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
  float r1, r2, r3;
  while(true)
  {
    /*
    This method ought to be cleaned up -- there is no reason for the while loop 
    since we're no longer accepting interactive user input.

    cout<<"\nPlease input 3 values for operational flexibility pairs(in nm) in increasing order:";
    cin>>r1>>r2>>r3;
    */
    r1 = ::atof(allInputs[currentInput++].c_str());
    r2 = ::atof(allInputs[currentInput++].c_str());
    r3 = ::atof(allInputs[currentInput++].c_str());
    if(!(r1>0 & r2>0 &r3>0 & r1<r2 & r2<r3)) { // valid values, all positive and in increasing order, then move to the next step
      cout<<"Operational flexbility values are invalid: Verify that r1 < r2 < r3.";
      exit(0);
    }
    else break;									// valid input
  }
  float* radii = new float[3];
  radii[0] = r1;
  radii[1] = r2;
  radii[2] = r3;
  routingDAG->generateOperFlexPairs(radii, 3, weatherData, deviationThreshold);	// generate the pairs of operational flexibility values	
  delete []radii;
  ctrl_OperFlexGenerated = OPER_FLEX_GENERATED;
  cout<<"\nOperational flexibility pairs successfully generated for the tree."<<endl;
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
    cout<<"\nPlease input the number of demands (entry nodes):";
    cin>>numDemands;
    if(numDemands>0)
      break;
    else
      cout<<"\nInvalid Value...";
  }
  string tempDemands;
  cin.ignore(1000 ,'\n');													// clear the buffer, ignore the leftover characters
  while(true)
  {
    cout<<"\nInput the demands one by one in nm (e.g. 2, 2, 1.5, 3, 2.7 means there are 5 entry nodes, and their rnp requirements are 2, 2, 1.5, 3 and 2.7.):";
    getline(cin, tempDemands);
    if(!inputDemandValid(tempDemands, numDemands))		// if the demand format is errorous
      cout<<"\nInvalid Input...";
    else if (!quadrant->demandFeasible(demandRNPs))		// the new demands are stored in the demandRNPs vector but cannot be accomodated by the quadrant
      cout<<"\nThe quadrant cannot accomodate the new demands...";
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
    cout<< endl << "Invalid weather count..." << endl;
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
  cout<<"\nReading weather data files now, please wait..."<<endl;
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
        cout<<"\nWeather not read in successfully..."<<endl;
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
      cout<<"\nWeather file directory error..."<<endl;
      return false;
    }
  } // loop over all weather files
  // weather files are read in, do some tests to make sure that the files are valid
  float totalProbabilityOfWeatherFiles = 0;
  for(int i=0; i<weatherData.size(); i++)
  {
    totalProbabilityOfWeatherFiles += weatherData[i].getProbability();
    cout << i << "   " << weatherData[i].getProbability() << endl;
  }
  // the total probability of the weather files is not 1
  if(abs(totalProbabilityOfWeatherFiles-1.0)>0.1)			
  {
    cerr<< "\nThe total probability of the weather data files is not 1."<<endl;
    cout << "Current Probability: " << totalProbabilityOfWeatherFiles << endl;
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
  //  cout << ctrl_QuadGenerated << endl << QUADRANT_GENERATED << endl;
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
  cout<< "\nTree Information successfully written to file.\n";
}
