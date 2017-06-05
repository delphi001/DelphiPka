/*
 * misc_timer.cpp
 *
 *  Created on: Apr 14, 2013
 *      Author: chuan
 */

#include "misc_timer.h"

#ifdef PARALLEL_OMP
#include <omp.h>  // LinWang : added for OpenMP timer 
#endif

using namespace std;


void CTimer :: formattedOutElapse(long int & tSec) const
{
   int tHour, tMin;

   if (tSec < 60) // less than 1 minute
   {
      
#ifdef PARALLEL_OMP
      double clockElapse = clockElapse_omp;  // LinWang : added for OpenMP timer 
#endif
 
      cout << fixed << (float)clockElapse/CLOCKS_PER_SEC << "   sec "<< endl;
   }
   else
   {
      tHour = tSec/3600; tSec = tSec - tHour*3600; 
      tMin  = tSec/60;   tSec = tSec - tMin*60;
   
      cout << tHour << ":" << tMin << ":" << tSec << endl;
   }

}


void CTimer :: start()
{   
   clockNow = clockStart;

#ifdef PARALLEL_OMP
   clockStart_omp = omp_get_wtime(); 
   clockNow_omp = clockStart_omp;   // LinWang : added for OpenMP timer 
#endif
   
   tNow = tStart; tmNowDateTime = tmStartDateTime;  
   cout << "\nProgram started on ";
   this->showTime(); 

}


void CTimer :: exit()
{    
   long int tSec;

   clockElapse = clock() - clockStart;

   tSec = clockElapse/CLOCKS_PER_SEC;

#ifdef PARALLEL_OMP

   clockNow_omp = omp_get_wtime();

   clockElapse_omp = (clockNow_omp - clockStart_omp)*CLOCKS_PER_SEC;

   tSec = clockElapse_omp/CLOCKS_PER_SEC;  // LinWang : added for OpenMP timer 

#endif
 
   tNow = time(0); tmNowDateTime = localtime(&tNow);  

   cout << endl << "   total CPU time ";

#ifdef PARALLEL_OMP
   cout << "with OpenMP on " << omp_get_num_procs() << " threads ";   // LinWang : added for OpenMP timer 
#endif

   cout << "is  ";

   this->formattedOutElapse(tSec); 

   cout << endl;

   cout << "Delphi exited on ";
   this->showTime();  

}


void CTimer :: showTime()
{
   tNow = time(0); tmNowDateTime = localtime(&tNow); 

   cout << 1900 + tmNowDateTime->tm_year << "-" << 1 + tmNowDateTime->tm_mon << "-" << tmNowDateTime->tm_mday << " at "
		<< tmNowDateTime->tm_hour << ":" << tmNowDateTime->tm_min << ":" << tmNowDateTime->tm_sec << endl;
}


void CTimer :: showElapse()
{   
   long int tSec;
      
   clock_t clockPast = clockNow;

   clockNow = clock(); 

   clockElapse = clockNow - clockPast;

   tSec = clockElapse/CLOCKS_PER_SEC;
   
#ifdef PARALLEL_OMP

   double clockPast_omp = clockNow_omp;

   clockNow_omp = omp_get_wtime();

   clockElapse_omp = (clockNow_omp - clockPast_omp)*CLOCKS_PER_SEC;  // LinWang : added for OpenMP timer 

   tSec = clockElapse_omp/CLOCKS_PER_SEC;

#endif
   
   this->formattedOutElapse(tSec);

}

