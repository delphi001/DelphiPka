#ifndef SPACE_EXCEPTIONS_H_
#define SPACE_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CZeroChargeRadius : public CWarning 
{
   public:
      CZeroChargeRadius(const int &cnt)
      {
         cerr << cnt << " CHARGED ATOMS HAVE ZERO RADIUS, WE CHANGED THEIR RADIUS TO 1" << endl;

      }
};


#endif /* SPACE_EXCEPTIONS_H_ */
