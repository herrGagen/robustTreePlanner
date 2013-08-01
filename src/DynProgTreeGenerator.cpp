#include "RoutingDAG.h"

#include "DynProgTreeGenerator.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>

/**
	\brief Anonymous namespace for exception handlers
*/
namespace {

class dynProgExceptions: public std::exception
{
public:
	static enum errorCodes
	{
		CHOOSE_TWO_ERROR = 1
	} errors;

	dynProgExceptions( unsigned int ERROR): error(ERROR){}
	unsigned int error;
	virtual const char *what() const throw()
	{
		switch(error)
		{
		case CHOOSE_TWO_ERROR:
			return "Bad argument(s) to chooseTwoIndex()";
		default:
			return "Unspecified exception in DynProgTreeGenerator.";
		}
	}
};

}

/**
	\brief ctor copying necessary info from a UserInterface instance.
*/
DynProgTreeGenerator::DynProgTreeGenerator(const UserInterface &UI)
{
	// First parse the demands table
	demands = UI.getDemandRNPs();

	// Second parse the number of nodes at each layer of the routingDAG
	// Store this in numInLayer, and its cumSum 
	RoutingDAG &graph = *(UI.getRoutingDAG() );
	numR = graph.getNumLayers();
	unsigned totalNodes = 0;
	for(unsigned int i = 0; i<graph.getNumLayers(); i++)
	{
		Node *layerChecker = graph.findNode(i,0);
		unsigned int j = 0;
		while(layerChecker != NULL)
		{
			layerChecker = graph.findNode(i,++j);
		}
		totalNodes += j;
		numInLayer.push_back(j);
		cumSumNumInLayer.push_back(totalNodes);
	}
	numT = numInLayer[0];
	resizeMemoTable( numT + totalNodes*(numT*numT - numT)/2 );

	// Next find all neighbor information and store in neighborsOf.
	neighborsOf.resize( totalNodes + 1 );
	double deviationThreshold = UI.getDeviationThreshold();
	double nodeEdgeThreshold  = UI.getNodeEdgeThreshold();
	const std::vector<WeatherData> &wDataSets = UI.getWeatherDataSets();
	double radius;
	// For each layer after the arrival node layer
	// Look at every node's neighbors and determine 
	// if the edge between them is safe.
	for(unsigned int i = 1; i<graph.getNumLayers(); i++)
	{		
		for(unsigned int j = 0; j<numInLayer[i]; j++)
		{
			Node *thisNode = graph.findNode(i,j);
			radius = thisNode->getDrawingRNP();
			unsigned int thisIndex = radiusAngleIndex(i,j);
			for(unsigned int k = 0; k< (unsigned int)thisNode->getNumInNeighbors(); k++)
			{				
				Edge *neighborEdge = thisNode->getInEdge(k);
				Node *neighborNode = thisNode->getInNode(k);
				bool nodeIsUnsafe = neighborEdge->isDangerousWeatherWithinLaneWidthW(radius, wDataSets, deviationThreshold, nodeEdgeThreshold);
				if( !nodeIsUnsafe )
				{
					unsigned int neighborIndex = neighborNode->getLayerIndex();
					neighborsOf[thisIndex].push_back(neighborIndex);
				}
			}		
		}
	}
	unsigned int sinkIndex = totalNodes;
	for(unsigned int fixIndex = 0; fixIndex < numInLayer[numR-1]; fixIndex++)
	{
		neighborsOf[sinkIndex].push_back(fixIndex);
	}
	generateTree();
}

/**
\brief Resizes the memoization table and fills with proper blank value.

Whenever we have not calculated the flow through a node, the memoization
table should contain the value -1 at the index associated with
that node / demand combination.

This function resizes the memoization table and properly fills it with
that value.
*/
void DynProgTreeGenerator::resizeMemoTable(unsigned int newSize)
{
	unsigned int oldSize = memoTable.size();
	if( newSize > oldSize )
	{
		memoTable.resize(newSize);
		std::fill(memoTable.begin()+oldSize, memoTable.end(), double(-1));
	}
}

/**
\brief Tells us if the flow(radInd,angInd,firstDemand,numDemands) has been calculated

\param radInd Radius index of node (think polar coordinates)
\param andIng Theta index of node
\param firstDemand Lowest angular index demand node we are covering
\param numDemands  Number of consecutive demand nodes we are covering
*/
bool DynProgTreeGenerator::isFlowMemoized(unsigned int radInd, 
	unsigned int angInd, 
	unsigned int firstDemand, 
	unsigned int numDemands) const
{
	if(numDemands == 0)
	{
		return true;
	}
	unsigned int index = getArrayIndex(radInd,angInd,firstDemand,numDemands);
	return( memoTable[index] > 0 );
}

/**
\brief Returns memoized flow(radInd,angInd,firstDemand,numDemands)

\param radInd Radius index of node (think polar coordinates)
\param andIng Theta index of node
\param firstDemand Lowest angular index demand node we are covering
\param numDemands  Number of consecutive demand nodes we are covering
*/
double DynProgTreeGenerator::getStoredFlow(unsigned int radInd, 
	unsigned int angInd, 
	unsigned int firstDemand, 
	unsigned int numDemands) const 
{
	if(numDemands == 0)
	{
		return 0;
	}
	unsigned int index = getArrayIndex(radInd,angInd,firstDemand,numDemands);
	return( memoTable[index] );

}

/**
\brief Turns 4 input ints into index into one dimensional vector

	Since alpha = (numT^2-numT)/2 is the total number of
	firstDemand,numDemands combinations

	And both radiusAngleIndex() and nChooseTwoIndex() run
	from 0 to their maximum values, a single-valued index into our 
	memoTable could be:

	radiusAngleIndex(radInd, angInd)*alpha + ...
	... nChooseTwoIndex(firstDemand, numDemands)

	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node
	\param firstDemand Lowest angular index demand node we are covering
	\param numDemands  Number of consecutive demand nodes we are covering
*/
unsigned int DynProgTreeGenerator::getArrayIndex(unsigned int radInd, 
	unsigned int angInd, 
	unsigned int firstDemand, 
	unsigned int numDemands) const
{

	// alpha in the discussion above
	unsigned int numCombos = (numT*numT - numT)/2;
	try
	{		
		unsigned int index = radiusAngleIndex(radInd, angInd)*numCombos + 
								nChooseTwoIndex(firstDemand,numDemands);
	return index;
	}
	catch( std::exception &e)
	{
		std::cerr << e.what() << "\nwith arguments: " << firstDemand << ", " << numDemands << std::endl;
		exit(-1);
	}	

}

/**
	\brief Translates pair (radInd, angInd) into a sequential index

	A number running from 0 to \sum_i numInLayer[i] - 1 is
	\sum_{i=0}^{radInd-1} numInLayer[i] + angInd

	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node
 
*/
unsigned int DynProgTreeGenerator::radiusAngleIndex(unsigned int radInd, unsigned int angInd) const
{

	unsigned int cumSum = radInd == 0 ? 0 : cumSumNumInLayer[radInd-1];
	unsigned int index = (cumSum + angInd);
	return index;

}

/**
	\brief Translates pair (firstDemand, numDemands) into a sequential index

	In our problem we have firstDemand in [0, numT-1]
	and numDemands in [1, numT-firstDemand]

	To get sequential indices for firstDemand (a) and numDemands (b) such that
	index( a, b ) + 1 = index(a, b + 1), we use the following simple formula:

	index = -1 + b + ... // simply making b the "ones digit"
	        + \sum_{i=0}^{a-1} (numT - i) // Number of feasible i,b combinations for i<a

	\param firstDemand Lowest angular index demand node we are covering
	\param numDemands  Number of consecutive demand nodes we are covering
	
*/
unsigned int DynProgTreeGenerator::nChooseTwoIndex( unsigned int firstDemand, unsigned int numDemands) const
{
	if(firstDemand >= numT ||
		numDemands > numT-firstDemand ||
		numDemands == 0 )
	{
		throw( dynProgExceptions( dynProgExceptions::CHOOSE_TWO_ERROR) );
	}
	unsigned int &i = firstDemand;
	unsigned int &j = numDemands;
	unsigned int index = numT*i - ( (i*i - i)/2 ) + j - 1;
	return index;
}

/**
	\brief Uses memoization table (or updates it) to return 

	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node
	\param firstDemand Lowest angular index demand node we are covering
	\param numDemands  Number of consecutive demand nodes we are covering

*/
double DynProgTreeGenerator::memoizedGetFlowThrough(unsigned int radInd, 
	unsigned int angInd, 
	unsigned int firstDemand, 
	unsigned int numDemands)
{

	if(numDemands ==0 )
	{
		return 0;
	}
	if(isFlowMemoized(radInd, angInd, firstDemand, numDemands) )
	{
		return getStoredFlow(radInd, angInd, firstDemand, numDemands);
	}
	else
	{
		return calculateAndStoreFlow(radInd, angInd, firstDemand, numDemands);
	}

}

/**
	\brief Calculates flow through node and stores result in memoization table.

	See class description for more detail, but implements
	f(r,t,a,b) = max[ max_{v' \in I(r,t)} f(v'.r, v'.t,a,b) ,
					  max_k (with {u' and w' in I(r,t)} ) of ...
					         ...  f(u.r, u.t, a, k) + f(w.r, w.t, a+k, b-k)

	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node
	\param firstDemand Lowest angular index demand node we are covering
	\param numDemands  Number of consecutive demand nodes we are covering
*/
double DynProgTreeGenerator::calculateAndStoreFlow(unsigned int radInd, 
	unsigned int angInd, 
	unsigned int firstDemand, 
	unsigned int numDemands)
{
	double highestFlow = 0;
	if(radInd == 0)
	{
		highestFlow = demands[angInd];
	}
	else
	{	
		// All incoming neighbors have radInd' = radInd-1, so we only need their angle indices
		const std::vector<unsigned int> &neighbors = getNeighbors(radInd,angInd);

		// max_{v' \in I(r,t)} f(v'.r, v'.t,a,b) 
		// i.e. find max flow from a single neighbor
		for( std::vector<unsigned int>::const_iterator nIter = neighbors.begin();
			nIter != neighbors.end(); 
			nIter++ )
		{
			double flow = memoizedGetFlowThrough(radInd-1, *nIter, firstDemand, numDemands);
			if(flow > highestFlow)
			{
				highestFlow = flow;
			}
		}

		// find max_k (with {u' and w' in I(r,t)} ) of ...
		//        ...  f(u.r, u.t, a, k) + f(w.r, w.t, a+k, b-k)
		// i.e. Find maximum flow achieved by merging two neighbor nodes at this node
		for(unsigned int k = 0; k< numDemands; k++)
		{
			for( std::vector<unsigned int>::const_iterator n1Iter = neighbors.begin();
				n1Iter != neighbors.end(); 
				n1Iter++ )
			{
				for( std::vector<unsigned int>::const_iterator n2Iter = n1Iter+1;
					n2Iter != neighbors.end(); 
					n2Iter++ )
				{
					double flow1 = memoizedGetFlowThrough(radInd-1, *n1Iter, firstDemand, k);
					double flow2 = memoizedGetFlowThrough(radInd-1, *n2Iter, firstDemand+ k, numDemands-k);
					if(flow1+flow2 > highestFlow)
					{
						highestFlow = flow1+flow2;
					}
				} // end for neighbor 1
			} // end for neighbor 2
		} // end for all values k at a split node
	} // end else not at an entry node
	memoizeFlow(highestFlow,radInd,angInd,firstDemand,numDemands);
	return highestFlow;

}

/**
	\brief Stores flow through node in memoization table.

	\param flow Value to be stored in table
	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node
	\param firstDemand Lowest angular index demand node we are covering
	\param numDemands  Number of consecutive demand nodes we are covering
*/
void DynProgTreeGenerator::memoizeFlow(double flow, unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands)
{
	unsigned int index = getArrayIndex(radInd,angInd,firstDemand,numDemands);
	memoTable[index] = flow;
}

/**
	\brief Returns all neighbors of the current node

	Currently all neighbors of the current node are all nodes with
	radInd' = radInd+1
	|angInd - angInd'| <= delta && angInd' \in [0, numT]

	\param radInd Radius index of node (think polar coordinates)
	\param andIng Theta index of node

	\retval angInd' of all neighbors, as the radInd' is always radInd+1
*/
const std::vector<unsigned int> &DynProgTreeGenerator::getNeighbors(unsigned int radInd, unsigned int angInd) const 
{
		
	unsigned int nodeIndex = radiusAngleIndex(radInd, angInd);
	return neighborsOf[nodeIndex];
}

/**
	\brief Fill memoization table
*/
void DynProgTreeGenerator::generateTree()
{

	for(unsigned int radInd = 0; radInd<numR; radInd++)
	{
		for(unsigned int angInd = 0; angInd < numInLayer[radInd]; angInd++)
		{
			std::cout << radiusAngleIndex(radInd, angInd) << '\t';
		}
	}
	std::cout << std::endl;

	memoizedGetFlowThrough(numR,0,0,numT);
#if 0
	for(unsigned int radInd = 0; radInd<numR; radInd++)
	{
		for(unsigned int angInd = 0; angInd < numInLayer[radInd]; angInd++)
		{
			for(unsigned int firstDemand = 0; firstDemand < numInLayer[radInd]; firstDemand++)
			{
				for(unsigned int numDemands = 1; numDemands < numInLayer[radInd]-firstDemand; numDemands++)
				{
					std::err << this->getArrayIndex(radInd,angInd,firstDemand,numDemands)
				}
			}
		}
	}
#endif
}