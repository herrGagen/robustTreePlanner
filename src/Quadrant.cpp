#include "Quadrant.h"
#include "math.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>

Quadrant::Quadrant(double cX, double cY, double ang, double iR, double oR, double iH, double oH)
{
	setcX(cX);
	setcY(cY);
	setAngle(ang);
	setiRadius(iR);
	setoRadius(oR);
	setiHeight(0);
	setoHeight(0);
	cHeight = oHeight - (oHeight-iHeight)*oRadius/(oRadius-iRadius);
	// inner variable, no setting function defined
	liftediRadius = sqrt(iRadius*iRadius+(iHeight-cHeight)*(iHeight-cHeight));
	liftedoRadius = sqrt(oRadius*oRadius+(oHeight-cHeight)*(oHeight-cHeight));
	liftStatus = QUADRANT_PLANE;
}

Quadrant::~Quadrant(void)
{
}

// modify and set the quadrant to a new configuration
void Quadrant::setQuadrant(double cX, double cY, double ang, double iR, double oR, double iH, double oH)
{
	setcX(cX);
	setcY(cY);
	setAngle(ang);
	setiRadius(iR);
	setoRadius(oR);
	setiHeight(iH);
	setoHeight(oH);
	cHeight = oHeight - (oHeight-iHeight)*oRadius/(oRadius-iRadius);
	// computed center height, inner variable of the class
	liftediRadius = sqrt(iRadius*iRadius+(iHeight-cHeight)*(iHeight-cHeight));
	liftedoRadius = sqrt(oRadius*oRadius+(oHeight-cHeight)*(oHeight-cHeight));
	if(iH!=0 || oH!=0)
		liftStatus = QUADRANT_LIFTED;
	else liftStatus = QUADRANT_PLANE; // both heights are 0, then its still on the plane z=0
}

// set the quadrant to its initial configuration
void Quadrant::reset()
{
	setQuadrant(0, 0, PI/6, 10, 35, 0, 0);
}

void Quadrant::setcX(double cX)
{
	centerX = cX;
}
void Quadrant::setcY(double cY)
{
	centerY = cY;
}

// normalize the angle to the range of [0, 2*PI)
void Quadrant::setAngle(double ang)
{
	if(ang<0)
	{
		int temp = 1;
		while(ang+2*PI*temp<0)
			temp++;
		angle = ang + 2*PI*temp;
	}
	else
	{
		int temp = 1;
		while(ang - 2*PI*temp > 0)
			temp++;
		angle = ang-2*PI*(temp-1);
	}
}

void Quadrant::setiRadius(double iR)
{
	if(iR>0)
		iRadius = iR;
	if(iR>oRadius)
		oRadius = iRadius+30;
}
void Quadrant::setoRadius(double oR)
{
	if(oR>0 && oR>iRadius)
		oRadius = oR;
	if(oR<iRadius)
		oRadius = iRadius+30;											// default size of the quadrant
}

void Quadrant::setiHeight(double iH)
{
	iHeight = iH;
	cHeight = oHeight - (oHeight-iHeight)*oRadius/(oRadius-iRadius);	// recompute the center height to help draw the quadrant
	liftediRadius = sqrt(iRadius*iRadius+(iHeight-cHeight)*(iHeight-cHeight));
	liftedoRadius = sqrt(oRadius*oRadius+(oHeight-cHeight)*(oHeight-cHeight));
}

void Quadrant::setoHeight(double oH)
{
	if(oH>=iHeight)
		oHeight = oH;
	else
		oHeight = iHeight+50;											// the default difference from iHeight to oHeight
	cHeight = oHeight - (oHeight-iHeight)*oRadius/(oRadius-iRadius);	// recompute the center height to help draw the quadrant
	liftediRadius = sqrt(iRadius*iRadius+(iHeight-cHeight)*(iHeight-cHeight));
	liftedoRadius = sqrt(oRadius*oRadius+(oHeight-cHeight)*(oHeight-cHeight));
}

double Quadrant::getcX()
{
	return centerX;
}
double Quadrant::getcY()
{
	return centerY;
}

double Quadrant::getAngle()
{
	return angle;
}

double Quadrant::getiRadius()
{
	return iRadius;
}

double Quadrant::getoRadius()
{
	return oRadius;
}

double Quadrant::getiHeight()
{
	return iHeight;
}

double Quadrant::getoHeight()
{
	return oHeight;
}

// set the lift status of the quadrant. If input status is an invalid value, do nothing
void Quadrant::setLiftStatus(int status)		
{
	if(status==QUADRANT_PLANE && iHeight == 0 && oHeight == 0)
		liftStatus = QUADRANT_PLANE;
	if(status==QUADRANT_LIFTED && (iHeight!=0 || oHeight!=0))							// at least one of them is not 0
		liftStatus = QUADRANT_LIFTED;
	liftediRadius = sqrt(iRadius*iRadius+(iHeight-cHeight)*(iHeight-cHeight));
	liftedoRadius = sqrt(oRadius*oRadius+(oHeight-cHeight)*(oHeight-cHeight));
}

int Quadrant::getLiftStatus()
{
	return liftStatus;
}

// decide if a group of demand for entry nodes is feasible or not, return true is feasilbe
bool Quadrant::demandFeasible(std::vector<double> &rnps)
{
	// first compute the maximum, minimum and sum of the rnp values	
	double rnpSum, minrnp, maxrnp;
	rnpSum = minrnp = maxrnp = rnps[0];
	for(unsigned int i=0; i<rnps.size(); i++)
	{
		rnpSum = rnpSum + 2*rnps[i];													// the sum of all the diameters
		minrnp = minrnp>rnps[i]? rnps[i] : minrnp;
		maxrnp = maxrnp<rnps[i]? rnps[i] : maxrnp;
	}
	if(rnpSum > PI*oRadius/2	||														// the outer arc length of the quadrant
		liftedoRadius - liftediRadius<=(log(double( rnps.size() ))/log(double(2))+1)*2*maxrnp ||	// there exists a way to merge the edges to a single fix node, in degree 2	 
		2*maxrnp> PI*iRadius/2)															// the inner boundary is long enough
		return false;
	return true;						// the demand is feasible
}

/*********************************************************************************************/
// functions used to generate routing graph structures
/*********************************************************************************************/
// the general function, which generated all the nodes in the routing DAG(directed graph)
// effectiveThres means if a weather cell's deviation probability is below this value, it's going to be considered as NULL
// routingThres means that we compute the weighted total probability of the weathercells p1*(0 or 1) + p2*(0 or 1) +..., 
// if the value < routingThres, then it is considered an obstacle 
bool Quadrant::generateDAG(std::vector<double> rnps, int n, double effectiveThres, double routingThres, const std::vector<WeatherData> &wDataSets, RoutingDAG* rDAG, double qAngleOffset, int numFixNodes)
{
	std::cout << "Generating DAG..." << std::endl;
	if(!generateEntryAndFixNodes(rnps, effectiveThres, routingThres, wDataSets, rDAG, qAngleOffset, numFixNodes))
		return false;
	std::cout << "Generating internal nodes..." << std::endl;

  generateRoutingDAGInternalNodes(rDAG, rnps, qAngleOffset);
	std::cout << "Finished generating internal nodes." << std::endl;
	return true;
}

// given a size n double array of entry points' rnps, and a weatherdata set, and the DAG structure that we are going to put the points in
bool Quadrant::generateEntryAndFixNodes(std::vector<double> rnps, double effectiveThres, double routingThres, const std::vector<WeatherData> &wDataSets, RoutingDAG *rDAG, double quadAngleOffset, unsigned int maxFixNodes)
{
  std::cout << "Max Number of Fix Nodes: " << maxFixNodes << std::endl;
	if(!demandFeasible(rnps))
	{
		std::cerr << "\nThe Quadrant is NOT large enough to be used for routing!"<<std::endl;
		return false;
	}
	/*******************************************************/
	// first generate a set of n entry nodes, equally distributed along the outer boundary of the quadrant
	double maxrnp = *max_element(rnps.begin(), rnps.end()); // the maximum rnp among the input rnps
	double startingAngle = 2*rnps[0]/oRadius; // the starting position and ending positions of the entry points
  double endingAngle = quadAngleOffset - 2*(*rnps.rbegin())/oRadius; // defined to be a little bit off the boundary of the quadrant
	for(unsigned int i=0; i<rnps.size(); i++)
	{
		double tempAngle = startingAngle + i*(endingAngle-startingAngle)/(rnps.size()-1);
		Node* tempNode = new Node(centerX+oRadius*cos(angle+tempAngle), centerY+oRadius*sin(angle+tempAngle), oHeight, ENTRY_NODE);
		tempNode->setLayer(0); 
		tempNode->setLayerIndex(i); // Entry Nodes are at layer 0 and layer index is set too
		rDAG->insertNode(tempNode); // insert the new Entry node into the routing graph
	}
	/***************************************************************************/
	// then generate a set of fix nodes, at least 1, at most 3
	double angleIncrement = (2*maxrnp+2)/iRadius;	
	// when we find a fix node, the next one should be at least 
	// angleIncrement aside from both sides

	double currentAngleLeft = quadAngleOffset/2+PI/18;									// iterate through boundary points on the left and right from the 45 degrees
  // WILLXYZ
	double currentAngleRight = quadAngleOffset/2;
	startingAngle = (maxrnp+0.5)/iRadius; // the boundary condition
	endingAngle = quadAngleOffset - (maxrnp+0.5)/iRadius;
	Node** nodeArrayToBeInsertedIntoTheDAG = new Node*[3]; 
	// store the fix nodes first, then insert them into the DAG in order of their angles

	double *anglesOfFixNodes = new double[3]; // store the angles of each fix nodes generated

	unsigned int numFixNodes = 0;
	while(currentAngleRight >= startingAngle || currentAngleLeft <= endingAngle)
	{
		if(currentAngleRight >= startingAngle)
		{
			Node* tempNode1 = new Node(centerX+iRadius*cos(angle+currentAngleRight), 
				centerY+iRadius*sin(angle+currentAngleRight), 
				iHeight, 
				FIX_NODE);

			// if not feasible
			if(tempNode1->isAnyWeatherCloserThanRadiusR(maxrnp, wDataSets, effectiveThres, routingThres))
			{
				delete tempNode1;
				currentAngleRight -= PI / 18; // test the next point on the inner arc ( 5 degrees out)
			}
			else // the node can serve as a fix node
			{
				nodeArrayToBeInsertedIntoTheDAG[numFixNodes] = tempNode1;
				// store the fix nodes information in the array

				anglesOfFixNodes[numFixNodes] = currentAngleRight;
				currentAngleLeft = std::max(currentAngleLeft, currentAngleRight + angleIncrement);
				currentAngleRight -= angleIncrement;
				numFixNodes++;
				if(numFixNodes>=maxFixNodes)
				{
					// already got enough fix nodes
					break;
				}
			}
		}
		if(currentAngleLeft<=endingAngle)
		{
			Node* tempNode2 = new Node(centerX+iRadius*cos(angle+currentAngleLeft), centerY+iRadius*sin(angle+currentAngleLeft), iHeight, FIX_NODE);
			if(tempNode2->isAnyWeatherCloserThanRadiusR(maxrnp, wDataSets, effectiveThres, routingThres))		// intersect with the weather
			{
				delete tempNode2;
				currentAngleLeft+=PI/18;
			}
			else																// the node can serve as a fix node
			{
				nodeArrayToBeInsertedIntoTheDAG[numFixNodes] = tempNode2;		// store the fix nodes information in the array
				anglesOfFixNodes[numFixNodes] = currentAngleLeft;
				currentAngleRight = std::min(currentAngleRight, currentAngleLeft - angleIncrement);
				currentAngleLeft += angleIncrement;
				numFixNodes ++;
				if(numFixNodes==3)												// already got enough fix nodes
				{
					break;
				}
			}
		}
	}
	 
		/* find at least 1 fix node,
		then insert them in order of the angle into the DAG,
		and set their layerIndex values */
	if(numFixNodes>=1)
	{
		/* we have the fix nodes stored in the nodeArrayToBeInsertedIntoTheDAG array, 
		and their corresponding angles are stores in the anglesOfFixNodes array */
		for(unsigned int i=0; i<numFixNodes; i++)	// insert the nodes one by one
		{
			int minAnglePos = 0;
			double minAngle = PI;
			for(unsigned int j=0; j<numFixNodes; j++)
			{
				if(anglesOfFixNodes[j]<minAngle)
				{
					minAngle = anglesOfFixNodes[j];
					minAnglePos = j; // find the position of the current minimum angle node
				}
			}
			anglesOfFixNodes[minAnglePos] = PI; // make it impossible to be the next minimum
			nodeArrayToBeInsertedIntoTheDAG[minAnglePos]->setLayerIndex(i);
			rDAG->insertNode(nodeArrayToBeInsertedIntoTheDAG[minAnglePos]); 
			// insert into the DAG based on increasing layerIndices
		}
		delete []nodeArrayToBeInsertedIntoTheDAG;
		delete []anglesOfFixNodes;
		return true;
	}
	delete []nodeArrayToBeInsertedIntoTheDAG;
	delete []anglesOfFixNodes;
	std::cerr << std::endl <<"The Quadrant is NOT large enough to be used for routing, please edit the quadrant!"<<std::endl;
	std::cout << std::endl << "Failure generating entry / fix nodes.!" << std::endl;
	return false;
}

// this function generated all the nodes in the DAG, but the edges are not generated, which have to be done by the generateEdgeSet() function defined in RoutingDAG.h
void Quadrant::generateRoutingDAGInternalNodes(RoutingDAG *rDAG, std::vector<double> rnps, double quadAngleOffset)
{
	unsigned int numLayers = 0;
	double minrnp = *max_element(rnps.begin(), rnps.end());
	std::cout << "minrnp initialized: " << minrnp << std::endl;
	for(unsigned int i=0; i<rnps.size(); i++)
	{
		if(rnps[i]!=0)
			minrnp = minrnp > rnps[i] ? rnps[i] : minrnp;							// use minrnp*3 as the interval in between levels, except 0
	}
	std::cout << "Found minrnp = " << minrnp << std::endl;
	if(liftedoRadius - liftediRadius > minrnp * 4)	
	{
		numLayers = 1;
		while((liftedoRadius - liftediRadius) / numLayers > minrnp * 3)
		{
			// compute the number of layers until the interval in between layer is less than maxrnp*3
			numLayers++;
		}
		numLayers--;
		// the actual number of layers was the number of divisions (which is 1 more)
	}
	std::cout << "Found numLayers = " << numLayers << std::endl;
	// generate the nodes on each layer, i is the layer index, where the smaller i denotes layer that is closer to the outer boundary
	for(unsigned int i=1; i<=numLayers; i++)
	{
		if (i % 5000 == 0) {
			std::cout << "Current layer is: " << i << std::endl;
		}
		// compute the unified heights and radius(on the base plane) of all the nodes on this layer
		double layerHeight = oHeight - (oHeight - iHeight) * i / (numLayers + 1);
		double layerRadius = oRadius - (oRadius - iRadius) * i / (numLayers + 1);

    // WILLXYZ It's possible that the following line needs to be updated to use quadAngleOffset
		double layerLength = PI*layerRadius/2;									// 1/4 circular perimeter
		int numNodesLayer = 0;													// numNodesLayer defines the number of nodes on the current layer
		while(numNodesLayer <= layerLength / (2 * minrnp + 2))
			numNodesLayer++;
		numNodesLayer--;														// the total number of nodes on current layer to ensure 2*maxrnp+2 interval
		double startingAngle = angle + quadAngleOffset / (2 * numNodesLayer);						// why 4? each angle is (PI/2)/numNodesLayer, but the starting and ending one are only half
		for(int j=0; j<numNodesLayer; j++)
		{
			double tempAngle = startingAngle+j*PI/(2*numNodesLayer);
			// compute the coordinates of the current node, set its layer and layerIndex information, then insert it into the DAG
			Node* tempNode = new Node(centerX+layerRadius*cos(tempAngle), centerY+layerRadius*sin(tempAngle), layerHeight, INTERNAL_NODE);
			tempNode->setLayer(i);
			tempNode->setLayerIndex(j);
			rDAG->insertNode(tempNode);
		}
	}
	rDAG->setNumLayers(numLayers+2);											// number of internal nodes layers + 2 layers for entry nodes and fix nodes
	rDAG->setNodesReadInStatus(NODES_READ_IN);
}

