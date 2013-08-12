#include "RoutingDAG.h"

#include "DynProgTreeGenerator.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>
#include <numeric>

/**
\brief Anonymous namespace for exception handlers
*/
namespace {
	/**
		\brief Exception class for this file

		Usage: throw(DynProgException( DynProgException::<errorCode> );

		try
		{
			codeThatThrowsError();
		}
		catch( std::exception &e )
		{
			std::cerr << e.what() << std::endl;
			... clean up code ..
		}
	*/
	class DynProgException: public std::exception
	{
	public:
		static enum errorCodes
		{
			RAD_INDEX_OUT_OF_BOUNDS,
			ANG_INDEX_OUT_OF_BOUNDS,
			FIRST_DEMAND_OUT_OF_BOUNDS,
			NUM_DEMANDS_OUT_OF_BOUNDS,
			NEIGHBORS_NOT_CALCULATED,
			NON_EXISTENT_NEIGHBOR
		} errors;

		DynProgException( errorCodes ERROR ): error(ERROR){}
		unsigned int error;
		virtual const char *what() const throw()
		{
			switch(error)
			{
			case RAD_INDEX_OUT_OF_BOUNDS:
				return "Radius index out of bounds";
			case ANG_INDEX_OUT_OF_BOUNDS:
				return "Angle index out of bounds";
			case FIRST_DEMAND_OUT_OF_BOUNDS:
				return "First demand out of bounds";
			case NUM_DEMANDS_OUT_OF_BOUNDS:
				return "Number of demands out of bounds";		
			case NEIGHBORS_NOT_CALCULATED:
				return "Neighbors not yet calculated";
			case NON_EXISTENT_NEIGHBOR:
				return "Requested neighbor does not exist";
			default:
				return "Unspecified exception in DynProgTreeGenerator.";
			}
		}
	};

}

/*****************************
Initialization Functions
*****************************/

/**
	\brief Simple ctor

	Makes all nodes neighbors of nodes in previous layers
	Makes all demand equal to 1

	\param numR Number of layers.
	\param numDem Number of demand nodes.
*/
DynProgTreeGenerator::DynProgTreeGenerator(unsigned int NUMR, unsigned int NUMDEM):
numR(NUMR),numDem(NUMDEM)
{
	// First set the demands vector.
	demands.resize(numDem);
	std::fill(demands.begin(),demands.end(),1);
	std::fill(numInLayer.begin(), numInLayer.end(), numDem);

	// Second parse the number of nodes at each layer of the routingDAG
	// Store this in numInLayer, and its cumSum 
	unsigned int totalNodes = 0;	

	for(unsigned int i = 0; i<numR; i++)
	{
		numInLayer.push_back(numDem);
		totalNodes += numDem;
		cumSumNumInLayer.push_back(totalNodes);
	}

	// Next find all neighbor information and store in neighborsOf.
	resizeMemoTable( numDem + totalNodes*(numDem*numDem + numDem)/2 );
	neighborsOf.resize(totalNodes + 1);

	// For each layer after the arrival node layer
	// Assume all are neighbors
	for(unsigned int i = 1; i<numR; i++)
	{		
		for(unsigned int j = 0; j<numDem; j++)
		{
			unsigned int thisIndex = radiusAngleIndex(i,j);
			for( unsigned int k = 0; k<numDem; k++ )
			{
				neighborsOf[thisIndex].push_back(k);
			}
		}
	}
	// Set sink node to be neighbor of all fix nodes
	numInLayer.push_back(1);
	cumSumNumInLayer.push_back(1+totalNodes);

	unsigned int sinkIndex = totalNodes;
	for(unsigned int fixIndex = 0; fixIndex < numInLayer[numR-1]; fixIndex++)
	{
		neighborsOf[sinkIndex].push_back(fixIndex);
	}
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
	numDem = numInLayer[0];
	resizeMemoTable( numDem + totalNodes*(numDem*numDem + numDem)/2 );

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
		// Loop over all nodes in current layer, 
		// Setting node(i,j) the active node
		for(unsigned int j = 0; j<numInLayer[i]; j++)
		{
			Node *thisNode = graph.findNode(i,j);
			radius = thisNode->getDrawingRNP();
			unsigned int thisIndex = radiusAngleIndex(i,j);
			// Find their neighbors and then 
			// determine if it's safe to go from neighbor[k] to the active node.
			for(unsigned int k = 0; k< (unsigned int)thisNode->getNumInNeighbors(); k++)
			{				
				Edge *neighborEdge = thisNode->getInEdge(k);
				Node *neighborNode = thisNode->getInNode(k);
				bool nodeIsUnsafe = neighborEdge->isDangerousWeatherWithinLaneWidthW(radius, wDataSets, deviationThreshold, nodeEdgeThreshold);
				if( !nodeIsUnsafe )
				{
					// This should match angInd for this neighbor.
					unsigned int neighborIndex = neighborNode->getLayerIndex();
					unsigned int neighborLayer = neighborNode->getLayer();
					if(i>0 && (neighborIndex >= numInLayer[i-1] || neighborLayer != i-1) )
					{
						continue; // THIS IS SKIPPING A BUG
					}
					neighborsOf[thisIndex].push_back(neighborIndex);
				}
			}		
		}
	}
	// Set sink node to be neighbor of all fix nodes
	numInLayer.push_back(1);
	cumSumNumInLayer.push_back(1+totalNodes);

	unsigned int sinkIndex = totalNodes;
	for(unsigned int fixIndex = 0; fixIndex < numInLayer[numR-1]; fixIndex++)
	{
		neighborsOf[sinkIndex].push_back(fixIndex);
	}
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
		usedNeighbors.resize(newSize);
		memoTable.resize(newSize);
		std::fill(memoTable.begin()+oldSize, memoTable.end(), double(-1));
	}
}

/*******************************
Memoization functions
*******************************/

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
	return( memoTable[index] >= 0 );
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

	if(numDemands == 0 )
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
	std::vector<unsigned int> predNeighbors;

	// start with an impossible flow to ensure memo entries with 0 flow properly
	// report their neighbors also have 0 flow
	double highestFlow = -1;
	int lowestAngleDifference = 10000;
	// If we are looking at a demand node, return a flow IFF
	// angle index matches starting index and number of demands is 1
	if(radInd == 0)
	{
		if(	angInd == firstDemand && numDemands == 1)
		{
			highestFlow = demands[angInd];
		}
		else
		{
			highestFlow = 0;
		}
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
			if(flow >= highestFlow ) 
			{
				// If tying an existing high flow, select neighbor with lowest turning angle
				if(flow == highestFlow)
				{
					int angleDifference = std::abs((int)*nIter - (int)angInd);
					if(angleDifference > lowestAngleDifference)
					{
						continue;
					}
					else
					{
						lowestAngleDifference = angleDifference;
					}
				}
				predNeighbors.resize(1);
				unsigned int neighborIndex = getArrayIndex(radInd-1,*nIter,firstDemand,numDemands);
				predNeighbors[0] = neighborIndex;
				highestFlow = flow;
			}
		}

		// Reset angle difference to prefer merging branches over single branches.
		lowestAngleDifference = 10000;
		// find max_k (with {u' and w' in I(r,t)} ) of ...
		//        ...  f(u.r, u.t, a, k) + f(w.r, w.t, a+k, b-k)
		// i.e. Find maximum flow achieved by merging two neighbor nodes at this node
		for(unsigned int k = 1; k< numDemands; k++)
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
					double flow2 = memoizedGetFlowThrough(radInd-1, *n2Iter, firstDemand + k, numDemands-k);
					if(flow1+flow2 >= highestFlow)
					{
						// If tying an existing high flow, select neighbor with lowest turning angle
						if(flow1+flow2 == highestFlow)
						{
							int angleDifference1 = (int)*n1Iter - (int)angInd;
							int angleDifference2 = (int)*n2Iter - (int)angInd;
							int angleDifference = std::abs(angleDifference1 + angleDifference2);
							if(angleDifference > lowestAngleDifference)
							{
								continue;
							}
							else
							{
								lowestAngleDifference = angleDifference;
							}
						}
						predNeighbors.resize(2);
						unsigned int n1Index = getArrayIndex(radInd-1, *n1Iter, firstDemand, k);
						predNeighbors[0] = n1Index;
						unsigned int n2Index = getArrayIndex(radInd-1, *n2Iter, firstDemand + k, numDemands-k);
						predNeighbors[1] = n2Index;						
						highestFlow = flow1+flow2;
					}
				} // end for neighbor 1
			} // end for neighbor 2
		} // end for all values k at a split node
	} // end else not at an entry node
	memoizeFlow(highestFlow,predNeighbors,radInd,angInd,firstDemand,numDemands);
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
void DynProgTreeGenerator::memoizeFlow(double flow, const std::vector< unsigned int> &predNeighbors, unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands)
{
	unsigned int index = getArrayIndex(radInd,angInd,firstDemand,numDemands);
	memoTable[index] = flow;
	usedNeighbors[index] = predNeighbors;
}

/*********************************
		Index Functions
**********************************/

/**
\brief Turns 4 input ints into index into one dimensional vector

Since alpha = (numDem^2-numDem)/2 is the total number of
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
	unsigned int numCombos = (numDem*numDem + numDem)/2;
	unsigned int index = radiusAngleIndex(radInd, angInd)*numCombos + 
		nChooseTwoIndex(firstDemand,numDemands);
	return index;
}

/**
	\brief returns radInd, angleInd from a tableIndex

	since tableIndex = radiusAngleIndex(radInd, angInd)*numCombos + 
					   nChooseTwoIndex(firstDemand,numDemands)

	And radiusAngleIndex(radInd, angInd) = cumSumNumInLayer[radInd-1] + angInd

	The inverse of this is:
	rI = radiusAngleIndex(radInd, angInd) = tableIndex / numCombos;
	radInd = index of first cumSum > rI
	angInd = rI - thisCumSum
*/
std::pair<unsigned int, unsigned int> DynProgTreeGenerator::radAngFromTableIndex( unsigned int tableIndex ) const
{
	unsigned int numCombos = (numDem*numDem + numDem)/2;
	unsigned int rI = tableIndex / numCombos;
	std::vector<unsigned int>::const_iterator rLayerIter; 
	rLayerIter = std::upper_bound(cumSumNumInLayer.begin(), cumSumNumInLayer.end(), rI );
	unsigned int radInd = rLayerIter - cumSumNumInLayer.begin();
	unsigned int angInd = rI;
	if(radInd > 0)
	{
			angInd = rI - *(rLayerIter - 1);
	}
	if(radInd > numR )
	{
		throw( DynProgException( DynProgException::RAD_INDEX_OUT_OF_BOUNDS) );
	}
	else if ( angInd >= numInLayer[radInd] )
	{
		throw( DynProgException( DynProgException::ANG_INDEX_OUT_OF_BOUNDS) );
	}
	return std::pair<unsigned int, unsigned int>(radInd, angInd);

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
	if(radInd > numR )
	{
		throw( DynProgException( DynProgException::RAD_INDEX_OUT_OF_BOUNDS) );
	}
	else if ( angInd >= numInLayer[radInd] )
	{
		throw( DynProgException( DynProgException::ANG_INDEX_OUT_OF_BOUNDS) );
	}
	unsigned int cumSum = radInd == 0 ? 0 : cumSumNumInLayer[radInd-1];
	unsigned int index = (cumSum + angInd);
	return index;

}

/**
\brief Translates pair (firstDemand, numDemands) into a sequential index

In our problem we have firstDemand in [0, numDem-1]
and numDemands in [1, numDem-firstDemand]

To get sequential indices for firstDemand (a) and numDemands (b) such that
index( a, b ) + 1 = index(a, b + 1), we use the following simple formula:

index = -1 + b + ... // simply making b the "ones digit"
+ \sum_{i=1}^{a} (numDem - i + 1) // Number of feasible i,b combinations for i<a

\param firstDemand Lowest angular index demand node we are covering
\param numDemands  Number of consecutive demand nodes we are covering

*/
unsigned int DynProgTreeGenerator::nChooseTwoIndex( unsigned int firstDemand, unsigned int numDemands) const
{
	if(firstDemand >= numDem )
	{
		throw( DynProgException( DynProgException::FIRST_DEMAND_OUT_OF_BOUNDS) );
	}
	else if( numDemands > numDem-firstDemand ||
		numDemands == 0 )
	{
		throw( DynProgException( DynProgException::NUM_DEMANDS_OUT_OF_BOUNDS) );
	}
	unsigned int &i = firstDemand;
	unsigned int &j = numDemands;
	unsigned int index = i*(numDem+1) - (i*i + i)/2 + j - 1;
	return index;
}

/**********************
Reporting Functions
**********************/

/**
\brief Returns all incoming neighbors of the current node

Currently all neighbors of the current node are all nodes with
radInd' = radInd - 1
|angInd - angInd'| <= delta && angInd' \in [0, numDem]

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
	\brief Fill memoization table iteratively (as opposed to recursively)
*/
void DynProgTreeGenerator::fillMemoTable()
{
	for(unsigned int r = 1; r<numR; r++)
	{
		for(unsigned int t = 0; t < numInLayer[r]; t++)
		{
			for(unsigned int a = 0; a < numDem; a++)
			{
				for(unsigned int n = 0; n<=numDem-a; n++)
				{
					memoizedGetFlowThrough(r,t,a,n);
				}
			}
		}
	}
	memoizedGetFlowThrough(numR,0,0,numDem);
}

/**
	\brief Recursively finds edges of tree used to reach this node in memoTable.

	Only call this manually after clearing the treeEdges vector.

	Recursion ends at the demand nodes which have no usedNeighbors

	\param tableIndex Root of the tree we will return
*/
void DynProgTreeGenerator::readTreeEdges( unsigned int tableIndex ) 
{

	const std::vector<unsigned int> &usedNeigh = usedNeighbors[tableIndex];

	for(std::vector<unsigned int>::const_iterator nIter = usedNeigh.begin();
		nIter != usedNeigh.end();
		++nIter)
	{
		treeEdges.push_back( std::pair<unsigned int, unsigned int>(tableIndex, *nIter) );
		readTreeEdges( *nIter);
	}
}

/**
\brief Saves edges in the best tree in the treeEdges vector
*/
void DynProgTreeGenerator::findBestTree( )
{
	double flow = getTotalFlow();
	treeEdges.clear();
	readTreeEdges( getArrayIndex(numR,0,0,numDem) );
}

/**
	\brief Copies best tree found to the DAG structure of the UI object

	\param UI The "UserInterface" object for the experiment.  Note this changes it.
*/
void DynProgTreeGenerator::writeBestTreeToDAG( UserInterface &UI )
{
	RoutingDAG *DAG = UI.getRoutingDAG();
	
	DAG->resetTree();
	for(unsigned int angInd = 0; angInd < numInLayer[0]; angInd++)
	{
		Node *thisNode = DAG->findNode(0, angInd);
		thisNode->setDrawingRNP(demands[angInd]);
	}
	findBestTree();
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator edgeIter;
	for( edgeIter = treeEdges.begin(); edgeIter != treeEdges.end(); ++edgeIter)
	{
		unsigned int highNodeIndex = edgeIter->first;
		unsigned int lowNodeIndex  = edgeIter->second;
		std::pair<unsigned int, unsigned int> highPair = radAngFromTableIndex(highNodeIndex);
		std::pair<unsigned int, unsigned int>  lowPair = radAngFromTableIndex( lowNodeIndex);
		Node *highRadNode = DAG->findNode(highPair.first, highPair.second );
		Node *lowRadNode  = DAG->findNode( lowPair.first,  lowPair.second );
		// The very first edge is the universal sink, something NOT in the DAG.
		// The DAG returns NUMM when it does not find a node.
		// Because of this, we skip all references to it in the tree.
		if( highRadNode == NULL )
		{
			continue;
		}
		Edge *theirEdge = DAG->findEdgeBetween(highRadNode, lowRadNode);
		highRadNode->setTreeNode();
		lowRadNode->setTreeNode();
		theirEdge->setTreeEdge();
		highRadNode->setDrawingRNP(lowRadNode->getDrawingRNP() );
		// the routingDAG expects nodes to report their used neighbor's position in their
		// personal neighbor list.  Not their actual index, so we have to do this bit of hackery.

		// Store the index of the neighbor used in the tree
		unsigned int neighborListIndex = 0;
		Node *highNodeFinder;
		for(unsigned int k = 0 ; k< lowRadNode->getOutSize(); k++)
		{
			highNodeFinder = lowRadNode->getOutNode(k);
			if(highNodeFinder == highRadNode)
			{
				neighborListIndex = k;
				break;
			}
		}
		if(highNodeFinder != highRadNode)
		{
			throw DynProgException( DynProgException::NON_EXISTENT_NEIGHBOR );
		}
		lowRadNode->setTreeOutEdgeIndex(neighborListIndex);
	}
	DAG->overrideTreeStatus();
	UI.inputOperationalFlexibility();
}

/**
\brief Prints edges in the tree

	Note: outputs current state, i.e findBestTree() should run first.

	Writes to cout using the format:
	(with node = destination node, and neig = source node)
	(NodeId, NeighborId)  node.radInd, node.angInd, neigh.radInd, neigh.angInd

*/
std::ostream & operator<<(std::ostream &os, const DynProgTreeGenerator& dp)
{
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator edgeIter;
	for( edgeIter = dp.treeEdges.begin(); edgeIter != dp.treeEdges.end(); ++edgeIter)
	{
		unsigned int me = edgeIter->first;
		unsigned int neighbor = edgeIter->second;
		os << "(" << me << ", " << neighbor << ")" << '\t';
		std::pair<unsigned int, unsigned int> myPair = dp.radAngFromTableIndex(me);
		std::pair<unsigned int, unsigned int> neighborPair = dp.radAngFromTableIndex(neighbor);
		os << myPair.first << ", " << myPair.second << ", ";
		os << neighborPair.first << ", " << neighborPair.second << std::endl;
	}
	return os;
}