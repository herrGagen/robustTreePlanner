#ifndef DYNPROGTREEGENERATOR_H
#define DYNPROGTREEGENERATOR_H

#include "UserInterface.h"
#include <vector>

/**
	\brief Creates a robust routing tree using dynamic programming

	Our layered DAG is indexed by integers r (radius) and t (for theta), with
	r=0 being the destination, and r=R_MAX being arrival nodes.

	We intend to calculate:
	f(r,t,a,b)

	v(r,t) node at (r,t)
	a Starting theta index of covered demand nodes
	b Number of consecutive demand nodes satisfied by node v(r,t)
	f() Number of arrivals in all demand nodes [a,a+b) that can be
	    safely routed through node v(r,t)

	This can be solved via Dynamic programming using the following recurrence:

	f(r,t,a,b) = max[ max_{v' \in I(r,t)} f(v'.r, v'.t,a,b) ,
					  max_k (with {u' and w' in I(r,t)} ) of ...
					         ...  f(u.r, u.t, a, k) + f(w.r, w.t, a+k, b-k)
	
	{I(r,t)} are the inbound neighbors of v(r,t)
*/
class DynProgTreeGenerator
{
public:	
	DynProgTreeGenerator(const UserInterface &ui);
	double memoizedGetFlowThrough(unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands);
	void generateTree();
protected:
	double calculateAndStoreFlow(unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands);
	double getStoredFlow(unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands) const;
	bool isFlowMemoized(unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands) const;
	void memoizeFlow(double flow, unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands);
	const std::vector<unsigned int> &getNeighbors(unsigned int radInd, unsigned int angInd) const;
private:
	unsigned int getArrayIndex(unsigned int radInd, unsigned int angInd, unsigned int firstDemand, unsigned int numDemands) const;
	unsigned int radiusAngleIndex(unsigned int radInd, unsigned int angInd) const;
	unsigned int nChooseTwoIndex( unsigned int firstDemand, unsigned int numDemands) const;
	void resizeMemoTable(unsigned int newSize);
private:
	/** This value used to test if we have calculated a flow value. */
	static const int NO_MEMO_ENTRY = -1;
	/** 
	\brief Memoization table holding maximum flow through nodes 

	Usage: This has to be sized at least |r|*|t|*(|a,b| < |t|^2/2)
	i.e. r t^3/2.

	Indices into this table are to be calculated using getArrayIndex.

	Note: For larger graphs, we will move to a sparse data representation 
*/
	std::vector<double> memoTable;
	unsigned int numR;   /**< Number of radius indices in layeredDAG */
	unsigned int numT;  /**< Total number of demand nodes in outer layer */
	std::vector<unsigned int> numInLayer;   /**< Number of nodes in each layer */
	std::vector<unsigned int> cumSumNumInLayer;   /**< cumSum of numInLayer */
	unsigned int delta;  /**< How many angle neighbors each node has */
	std::vector< std::vector< unsigned int > > neighborsOf;
	std::vector<double> demands;
};

#endif