/*
 * energy_exceptions.h
 *
 *      Author: Lin Wang, lwang3@g.clemson.edu
 *
 * These warnings and exceptions should inherit from the classes CWarning and CException in order to
 * keep track the number of warnings and maintain consistent format.
 */

#ifndef ENERGY_EXCEPTIONS_H_
#define ENERGY_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CIonicCalcUnderTest : public CWarning
{
   public:
	CIonicCalcUnderTest()
      {
         cerr << "THIS PART IS STILL UNDER TESTING. \n";
      }
};

class CReacFieldEnergyOverEST : public CWarning
{
   public:
	CReacFieldEnergyOverEST()
      {
         cerr << "BE CAREFUL!! REACTION FIELD ENERGY MIGHT BE OVERESTIMATED!!!! \n";
      }
};

class CThreadsLessThanProcs : public CWarning
{
	public:
	 CThreadsLessThanProcs(const int& iProcs, const int&iThreads)
		{
			cerr << "PROGRAM DETECTS " << iProcs << " PROCESSORS AVAILABLE. HOWEVER, PROGRAM IS SET " << iThreads << " THREADS FOR CALCULATION. \n";
			cerr << "FOR BEST PERFORMANCE, SET THREADS NUMBER NO LESS THAN TOTAL NUMBER OF PROCESSORS BY USING COMMAND 'export OMP_NUM_THREADS=Processors_Number'. \n";
		}	
};

#endif /* ENERGY_EXCEPTIONS_H_ */
