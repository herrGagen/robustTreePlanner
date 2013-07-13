#include "DemandProfile.h"
#include <fstream>
#include "math.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <typeinfo>

DemandProfile::DemandProfile(void)
{
  numDemands = 0;
  minX = maxX = minY = maxY = 0;
}

DemandProfile::~DemandProfile(void)
{
  xCoors.clear();
  yCoors.clear();
  rnp.clear();
}

// read in the demand profile from file, the content is passed in by the character array buffer, and the size is fileSize
bool DemandProfile::readInFile(std::string fileName)			
{

  std::cout << "Parsing demand file: " << fileName << std::endl;
  ifstream dataStream( fileName.c_str(), ifstream::in );
  dataStream.seekg(0, ios_base::end);
  dataStream.seekg(0, ios_base::beg);
  std::string thisLine;

  // the first step tries to read in the center location(airport location)
  rnp.clear();
  xCoors.clear();
  yCoors.clear();
  xCoors.reserve(128);
  yCoors.reserve(128);

  /*
    The data files are formatted in the following way (# start comment lines):	
    Then file has the following parameters:

    CASE_NAME ==> Ignore
    CURRENT_TIME ==> Ignore
    SIMULATION_TIME_RANGE ==> timeStart [SPACE] timeEnd
    CENTER_LOCATION ==> centerX [SPACE] centerY
    AIRPORT_RANGE_RING_IN_NAUTICAL_MILES ==> Ignore
    RNP_LEVELS ==> Ignore
    RNP_DISTRIBUTION ==> Ignore

    ### Flight Data in next comment block
  */
  while ( !dataStream.eof() )
    {
      size_t found = 0;
      if( dataStream.peek() == '#' || dataStream.peek() == ' ')
        {
          std::getline(dataStream, thisLine);
          continue;
        }
      std::getline(dataStream, thisLine);

      // Look for SIMULATION_TIME_RANGE and store following 2 numbers in timeStart, timeEnd
      if( (found = thisLine.find("SIMULATION_TIME_RANGE") ) != std::string::npos )
        {
          size_t beginNumber = thisLine.find_first_of("-.0123456789",found);
          size_t endNumber = thisLine.find_first_not_of("-.0123456789", beginNumber);
	
          timeStart = (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str(); // store it in the timeStart std::string
          beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
          endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
          endNumber = (endNumber < thisLine.length() ) ? endNumber : thisLine.length()-1;
          timeEnd   = (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str(); // store it in the timeStart std::string
        }

      // Look for CENTER_LOCATION and store following 2 numbers in centerX, centerY
      if( (found = thisLine.find("CENTER_LOCATION")) != std::string::npos )	
        {
          size_t beginNumber = thisLine.find_first_of("-.0123456789",0);
          size_t endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
          centerX = (double) ::atof( (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str() ); // store it in the timeStart std::string
          beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
          endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
          endNumber = (endNumber < thisLine.length() ) ? endNumber : thisLine.length()-1;
          centerY   = (double) ::atof( (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str() ); // store it in the timeStart std::string
        }

      // Look for BEGIN_FLIGHTS to signal that we're in the correct location for flight data
      if( (found = thisLine.find("BEGIN_FLIGHTS") ) != std::string::npos )	
        {
          break;
        }

    }

  /*
    Flight data is in the next portion of the file, between the "tags" BEGIN_FLIGHTS, END_FLIGHTS
    and is in the format:
    Unique ID,ACID,RNP Level,Airport Range Entry Time,Airport Range Entry Lat,Airport Range Entry Lon, Dept. Airport,Arr. Airport,Flight Plan std::string
    Unique ID ==> Ignore
    ACID ==> Ignore
    RNP Level ==> std::vector<double>  rnp (note, the entries are all integers)
    Airport Range Entry Time ==> Ignore
    Airport Range Entry Lat ==> std::vector xCoors
    Airport Range Entry Lon ==> std::vector yCoors
    Dept. Airport ==> Ignore
    Arr. Airport ==> Ignore
    Flight Plan std::string ==> Ignore
  */
  std::string tempString;
  double tempdouble = 0;
  int numPaths = 0;
  while ( !dataStream.eof() )
    {
      // Ignore comment lines
      if( dataStream.peek() == '#')
        {
          std::getline(dataStream, thisLine);
          continue;
        }
      std::getline(dataStream, thisLine);

      // Break after we're done with flight data
      if( thisLine.find("END_FLIGHTS") != std::string::npos )	
        {
          break;
        }

      // Ignore the first 2 comma separated values
      unsigned long beginNumber = thisLine.find_first_of(",",0); // find first comma
      unsigned long endNumber = thisLine.find_first_of(",",beginNumber+1); // find second comma

      // Find the RNP value
      beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
      endNumber   = thisLine.find_first_not_of("-.0123456789",beginNumber);
		
      tempString = (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str();
      tempdouble = (double) ::atof( tempString.c_str() );
      rnp.push_back( tempdouble );

      // Now ignore next CSV
      beginNumber = thisLine.find_first_of(",",endNumber); // first comma after RNP
      endNumber   = thisLine.find_first_of(",",beginNumber+1); // find second comma after RNP

      // The next two numbers are the X and Y coordinates
      // First the xCoord
      beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
      endNumber   = thisLine.find_first_not_of("-.0123456789",beginNumber);
					
      tempString = (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str();
      tempdouble = (double) ::atof( tempString.c_str() );
      xCoors.push_back( (double)tempdouble );

      // Next we find the yCoord
      beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
      endNumber   = thisLine.find_first_not_of("-.0123456789",beginNumber);
      endNumber = endNumber < thisLine.length() ? endNumber : thisLine.length();

      tempString = (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str();
      tempdouble = (double) ::atof( tempString.c_str() );
      yCoors.push_back( (double)tempdouble );
    }
  return handleInputData();
}


bool DemandProfile::handleInputData()
{
  if(centerX>360 || centerY>360)							// the lati/long of the center is invalid
    {
      return false;
    }
  if(centerX>180)	
    {
      centerX -=360;							// normalize the lati/longs to (-180, 180] range
    }
  if(centerY>180)	
    {
      centerY -=360;
    }
  numDemands = xCoors.size();								// the total number of demands(flights) in the profile
  // the demand profile is almost on a circle centers at position(centerX, centerY)
  minX = maxX = centerX;
  minY = maxY = centerY;
  // check if there is input content error, and regulate the lati/long to the range of (-180, 180]
  for(unsigned int i=0; i<numDemands; i++) 
    {
      if(rnp[i]<0 || xCoors[i]>360 || yCoors[i]>360)
        {
          std::cout << rnp[i] << " " << xCoors[i] << " " << yCoors[i] << std::endl;
          std::cerr << "\nFile Content Error!"<<std::endl;				// prompt that the file has format errors
          return false;
        }
      if(xCoors[i]>180)									// change the format to (-180, 180]
        {
          xCoors[i]-=360;
        }
      if(yCoors[i]>180)
        {
          yCoors[i]-=360;
        }
      if(xCoors[i]<minX)	
        {
          minX = xCoors[i];
        }
      if(yCoors[i]<minY)	
        {
          minY = yCoors[i];
        }
      if(xCoors[i]>maxX)	
        {
          maxX = xCoors[i];
        }
      if(yCoors[i]>maxY)	
        {
          maxY = yCoors[i];
        }
    }
  // find the latitude diffference
	// the demand circle has total diameter of 200 pixels(a given value, best value chosen for display purpose) 
	// Joe changed the next 2 lines, March 6, 2013, units are degrees per pixel, assuming 200nm outer radius, 3438nm radius of earth
	latiPerPixel =  (200*(180/PI)/3438)/TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL;   //   (maxX - minX)/(2*TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL);						
	longPerPixel =  (200*(180/PI)/(3438*cos(centerX*(PI/180))))/TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL;   //   (maxY - minY)/(2*TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL);
  return true;
}

void DemandProfile::reset()
{
  xCoors.clear();
  yCoors.clear();
  rnp.clear();
  numDemands = 0;
  minX = maxX = minY = maxY = 0;
}


// thsi function is called when reading in weather data. The function is used to do data conversion, from
// the read in weather cell in lati/long to their x/y positions. 
void DemandProfile::getDemandInfo(double* cX, double* cY, double* laPerPixel, double* loPerPixel, std::string* timeS, std::string* timeE)
{
  *cX = centerX;
  *cY = centerY;
  *laPerPixel = latiPerPixel;
  *loPerPixel = longPerPixel;
  *timeS = timeStart;
  *timeE = timeEnd;
}

// The function is also used to do weather data triming,
// where the data cells that are way too far from the demand are not stored in the memory at all
void DemandProfile::getRange(double *minLati, double *minLong, double *maxLati, double *maxLong)
{
  *minLati = minX;	*minLong = minY;
  *maxLati = maxX;	*maxLong = maxY;
}

// generate a set of 8 entry points and their rnp values based on the starting angle of the quadrant
void DemandProfile::generateDemandVector(std::vector<double> &demand, double startingAngle, double endingAngle, double nmilesPerPixel)
{
  double demandRNPS[NUM_ENTRY_NODES_PER_QUADRANT] = {0};
  double angleIncrement = (endingAngle - startingAngle)/NUM_ENTRY_NODES_PER_QUADRANT;
  for(unsigned int i=0; i<numDemands; i++)
    {
      // compute the x, y coordinates in the opengl coordinate system for each demand, where the center is (0, 0)
      double tempX = (xCoors[i]-centerX)/latiPerPixel;	
      double tempY = (yCoors[i]-centerY)/longPerPixel;
      if(tempX!=0)
        {
          double angle = atan(tempY/tempX);					// the center angle of the point
          if(tempX<0)	//if the point is in the 2nd or the 3rd quadrant
		  {
            angle = angle+PI;
		  }
          else if(tempX>0 && tempY<0) // in the 4th quadrant
		  {
            angle = angle+2*PI;	 // make it in the 0--2*PI range
		  }
          // now we have the angle in range (0, 2*PI)
          if(endingAngle<2*PI)	
            {
              if(angle>startingAngle && angle<=endingAngle)
                {
                  int tempIndex = int((angle-startingAngle)/angleIncrement);
                  demandRNPS[tempIndex] = std::max(demandRNPS[tempIndex], rnp[i]);	 // compute the std::max rnp value for this range of entry nodes
                }
            }
          else	// starting and ending angle are on different sides of the positive x axis
            {
              if(angle>startingAngle || angle<endingAngle-2*PI)
                {
                  int tempIndex = angle>=startingAngle ? int((angle-startingAngle)/angleIncrement) :  int((angle+2*PI-startingAngle)/angleIncrement);
                  demandRNPS[tempIndex] = std::max(demandRNPS[tempIndex], rnp[i]);	 // compute the std::max rnp value for this range of entry node
                }
            }			
        }
    }
  // convert to opengl coordinate system and in the project, rnp is half the lane width
  for(unsigned int i=0; i<NUM_ENTRY_NODES_PER_QUADRANT; i++)		
    {
    demand.push_back(demandRNPS[i]/(2*nmilesPerPixel));						// store them into the resulting std::vector
    }
}
