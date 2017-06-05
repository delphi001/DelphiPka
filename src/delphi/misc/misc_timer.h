/**
 * @file misc_timer.h
 * @brief a timer class to print program execution time
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <iostream>
#include <time.h>  /* or include ctime */
#include "../interface/interface_abstractmodule.h"

class CTimer
{
   private:
      //for print out current date and time
      time_t tStart, tNow;
      tm * tmStartDateTime, * tmNowDateTime;
      //for precise calculating elapsed time
      clock_t clockStart, clockNow, clockElapse; 

	  double clockStart_omp, clockNow_omp, clockElapse_omp; // LinWang : added for OpenMP timer 
      
      void formattedOutElapse(long int &) const;

   public:

      /**
       * constructor
       */
      CTimer()
      { 
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                   CTimer is constructed                      *\n";
         cout << "****************************************************************\n";
#endif

         tStart = time(0); tmStartDateTime = localtime(&tStart); 
         clockStart = clock();

      };
      
      /**
       * destructor
       */
      ~CTimer()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                    CTimer is destroyed                       *\n";
         cout << "****************************************************************\n";
#endif
      };
      
      /**
       * function to output starting time of the program
       */
      void start();

      /**
       * function to output exiting time of the program
       */
      void exit();
      
      /**
       * function to show current time
       */
      void showTime();
      
      /**
       * function to show elapsed time from last call
       */
      void showElapse(); 
};


#endif // TIMER_H_
