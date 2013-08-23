#include <iostream>

#include "UserInterface.h"
#include "TreeVerifier.h"
#include "DynProgTreeGenerator.h"

int main(int argc, char* argv[])
{	

  std::string inputFile = "inputs.txt";
  if (argc > 1)
  {
    for (int i = 0; i < argc; i++)
    {
      std::string argvi = argv[i];
      std::cout << argvi << std::endl;
      if (argvi == "-cinput")
      {
        inputFile = argv[i+1];
      }
    }
  }
  
	UserInterface RTP;
	if(RTP.ProgramBegins(inputFile) )
	{
		RTP.makeRTPTreeAndFinish();
		TreeVerifier checker(RTP);
		checker.appendReportToFiles("InvalidTrees.txt","ValidTrees.txt");
	    checker.outputTreeLengthStats();
	}

	UserInterface dpDagMaker;
	if( dpDagMaker.ProgramBegins(inputFile) )
	{
		DynProgTreeGenerator dp(dpDagMaker);
		dp.writeBestTreeToDAG(dpDagMaker);
		TreeVerifier dpChecker(dpDagMaker);
		dpChecker.appendReportToFiles("dpInvalidTrees.txt","dpValidTrees.txt");
		dpChecker.outputTreeLengthStats();
		dpDagMaker.setOutputFileName( dpDagMaker.getOutputFileName() + ".dp.xml" );
		dpDagMaker.saveTreeInformation();
	}

	return 0;
}
