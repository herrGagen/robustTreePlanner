#include "DynProgTreeGenerator.h"
#include "UserInterface.h"
#include <exception>
#include <memory>
#include <algorithm>
#include <numeric>

#include "DynProgTreeGeneratorTester.h"

DynProgTreeGeneratorTester::DynProgTreeGeneratorTester(  std::string inputFileName )
{

	UserInterface userInterface;
	userInterface.ProgramBegins(inputFileName);
	generator = DynProgTreeGenerator(userInterface);
}

DynProgTreeGeneratorTester::DynProgTreeGeneratorTester( unsigned int NUMR, unsigned int NUMDEM):
generator( DynProgTreeGenerator(NUMR,NUMDEM) )
{
}

/**
	\brief Ensure that array index errors if 0 is passed as number of demands
*/
bool DynProgTreeGeneratorTester::indexFailsOnZeroDemands() const
{

	bool failsOnZero = false;
	try 
	{
		generator.getArrayIndex(0,0,0,0);
	}
	catch( ... )
	{
		failsOnZero = true;
	}
	return failsOnZero;
}

/**
	\brief Ensure that all indices are sequential.
*/
bool DynProgTreeGeneratorTester::indexCountsSequentially() const 
{
	unsigned int numR = generator.numR;
	unsigned int numDem = generator.numDem;
	bool retval = true;
	unsigned int lastIndex = 0;
	for(unsigned int i = 0; i<numR; i++)
	{		
		for(unsigned int j = 0; j<generator.numInLayer[i]; j++)
		{
			for(unsigned int a = 0; a<numDem; a++)
			{
				for(unsigned int b = 1; b<=numDem-a; b++)
				{
					unsigned int index = generator.getArrayIndex(i,j,a,b);
#if defined(INDICES_NOT_WORKING)
					std::cout << index << '\t';
#endif
					if(index > 0)
					{
#if !defined(INDICES_NOT_WORKING)
						EXPECT_EQ(index - lastIndex, 1) << index << ", " << lastIndex << '\n' << i << ", " << j << ", " << a << ", " << b << std::endl;
#endif
						if(index - lastIndex != 1)
						{
							retval = false;
						}
					}
					lastIndex = index;
				}
#if defined(INDICES_NOT_WORKING)
				std::cout << "( " << i << ", " << j << ", " << a << ")" << std::endl;
#endif
			}
		}
	}	
	return retval;
}

/**
	\brief Checks inverse indexing function
*/
bool DynProgTreeGeneratorTester::reverseIndexVerifier() const 
{
	unsigned int numR = generator.numR;
	unsigned int numDem = generator.numDem;
	bool retval = true;
	for(unsigned int i = 0; i<numR; i++)
	{		
		for(unsigned int j = 0; j<generator.numInLayer[i]; j++)
		{
			for(unsigned int a = 0; a<numDem; a++)
			{
				for(unsigned int b = 1; b<=numDem-a; b++)
				{
					unsigned int index = generator.getArrayIndex(i,j,a,b);
					std::pair<unsigned int, unsigned int> radAng = generator.radAngFromTableIndex( index );
					EXPECT_EQ(i,radAng.first)  << "radIndex inverse equality check (a,b) (" << a << ", " << b << ")";
					EXPECT_EQ(j,radAng.second) << "angIndex inverse equality check (a,b) (" << a << ", " << b << ")";
					if(radAng.first != i || radAng.second != j)
					{
						retval = false;
					}
				}
			}
		}
	}	
	return retval;
}

std::pair<double,double> DynProgTreeGeneratorTester::getInputAndOutput() 
{
	const std::vector<double> &dem = generator.demands;
	double allDemand = std::accumulate(dem.begin(), dem.end(), (double)0.00 );
	double dpSatisfiedDemand = generator.getTotalFlow();
	return std::pair<double,double>(allDemand, dpSatisfiedDemand );
}

unsigned int DynProgTreeGeneratorTester::numBranchesNotTerminatingAtArrivalNodes()
{
	unsigned int retval = 0;
	generator.findBestTree();
	const std::vector<std::pair<unsigned int, unsigned int> > &treeEdges = generator.treeEdges;
	std::vector<std::pair<unsigned int, unsigned int> >::const_iterator edgeIter;
	for( edgeIter = treeEdges.begin()+1; edgeIter != treeEdges.end(); ++edgeIter)
	{
		if(edgeIter->first > (edgeIter-1)->first )
		{
			std::pair<unsigned int, unsigned int> lastRA = generator.radAngFromTableIndex((edgeIter-1)->second);
			std::pair<unsigned int, unsigned int> thisRA = generator.radAngFromTableIndex(edgeIter->second);

			if( lastRA.first != 0 )
			{
				retval++;
			}
		}
	}
	return retval;
}


////////// Implement actual tests ///////////////////////

TEST_F(DynProgTreeGeneratorTester_TXT, DemandConservation)
{

	std::pair<double, double> retVals = getInputAndOutput();
	ASSERT_PRED_FORMAT2(::testing::DoubleLE, retVals.second,  retVals.first ) << "DP output does not match input demand.";
	EXPECT_DOUBLE_EQ( retVals.first, retVals.second) << "Nonfatal: DP does not allow all flow.";
	unsigned int badBranchCount = numBranchesNotTerminatingAtArrivalNodes();
	EXPECT_EQ( badBranchCount, 0 ) << "Bad branches found.";
}

TEST_F(DynProgTreeGeneratorTester_TXT, IndexTests)
{
	EXPECT_TRUE( indexCountsSequentially() ) << "Non-consecutive index found.";
	EXPECT_TRUE( indexFailsOnZeroDemands() ) << "0 NumDemands not returning 0.";
	EXPECT_TRUE( reverseIndexVerifier() ) << "Reversing index encoding not working";
}

TEST_F(DynProgTreeGeneratorTester_Default, DemandConservation)
{
	std::pair<double, double> retVals = getInputAndOutput();
	ASSERT_PRED_FORMAT2(::testing::DoubleLE, retVals.second,  retVals.first) << "DP output greater than input demand.";
	EXPECT_DOUBLE_EQ( retVals.first, retVals.second) << "Nonfatal: DP does not route all traffic.";
	unsigned int badBranchCount = numBranchesNotTerminatingAtArrivalNodes();
	EXPECT_EQ( badBranchCount, 0 ) << "Bad branches found.";
}

TEST_F(DynProgTreeGeneratorTester_Default, IndexTests)
{
	EXPECT_TRUE( indexFailsOnZeroDemands() ) << "0 NumDemands not returning 0.";
	EXPECT_TRUE( indexCountsSequentially() ) << "Non-consecutive index found.";
	EXPECT_TRUE( reverseIndexVerifier() ) << "Reversing index encoding not working";
}


