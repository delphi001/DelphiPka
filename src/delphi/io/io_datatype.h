/**
 * @file io_datatype.h
 * @brief defining data types used for IO class
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef IO_DATATYPE_H_
#define IO_DATATYPE_H_

#include <string.h>

#include "../interface/environment.h"
#include "../misc/misc_grid.h"
#include "io_exceptions.h"

/**
 * Class to store force field read from .siz and .crg files. Notice that values are all defined to be
 * private so that accessing these values have to be via the provided public member functions.
 */
class CForce
{
   private:
      string strAtom;       // atom name          a6
      string strResidue;    // residue name       a3
      string strResidueNum; // residue number     a4
      string strChain;      // subunit name       a1
      delphi_real fValue;   // atom radii/charge  f8.4
     
   public:
      /**
       * default constructor
       */
      CForce()
      {
         this->strAtom       = " ";
         this->strResidue    = " ";
         this->strResidueNum = " ";
         this->strChain      = " ";
         this->fValue        = 0.0;
      };

      /**
       * function to set atom name
       * @param strAtomIn atom name read from input file
       */
      void setAtom(const string strAtomIn) {this->strAtom = strAtomIn;}
      
      /**
       * function to set residue name
       * @param strResidueNameIn residue name read from input file
       */
      void setResidue(const string strResidueNameIn) {this->strResidue = strResidueNameIn;}
      
      /**
       * function to set residue number
       * @param strResidueNumIn residue number read from input file
       */
      void setResidueNum(const string strResidueNumIn) {this->strResidueNum = strResidueNumIn;}
  
      /**
       * function to set chain name
       * @param strChainIn chain name read from input file
       */
      void setChain(const string strChainIn) {this->strChain = strChainIn;} 
      
      /**
       * function to set the value (either charge or radius) associated with a particular atom
       * @param fValueIn force field value associated with this atom
       */
      void setValue(const delphi_real fValueIn) {this->fValue = fValueIn;}
  
      /**
       * function to get atom name
       * @return atom name
       */
      string getAtom() const {return this->strAtom;}

      /**
       * function to get residue name
       * @return residue name
       */
      string getResidue() const {return this->strResidue;}
      
      /**
       * function to get residue number
       * @return residue number
       */
      string getResidueNum() const {return this->strResidueNum;}
      
      /**
       * function to get chain name
       * @return chain name
       */
      string getChain() const {return this->strChain;}
      
      /**
       * function to get force field value
       * @return charge or radius of this atom
       */
      delphi_real getValue() const {return this->fValue;}
};

/**
 * class to store atom info read from .pdb file. Notice that values are all defined to be
 * private so that accessing these values have to be via the provided public member functions.
 */
class CAtomPdb // delphi_pdb_file_record
{
   private:
      delphi_real fRadius;
      delphi_real fCharge;
      SGrid<delphi_real> gPose;
      string strAtInf;

   public:
      /**
       * default constructor
       */
      CAtomPdb()
      {
         this->fRadius  = 0.0;
         this->fCharge  = 0.0;
         this->gPose.nX = 0.0; this->gPose.nY = 0.0; this->gPose.nZ = 0.0;
         this->strAtInf = " ";
      }; 

      /**
       * function to set radius of this atom
       * @param fRadiusIn radius of this atom
       */
      void setRadius(const delphi_real fRadiusIn) {this->fRadius = fRadiusIn;}
      
      /**
       * function to set charge of this atom
       * @param fChargeIn charge on this atom
       */
      void setCharge(const delphi_real fChargeIn) {this->fCharge = fChargeIn;}
      
      /**
       * function to set position of this atom
       * @param gPoseIn position of this atom
       */
      void setPose(const SGrid<delphi_real> gPoseIn) {this->gPose = gPoseIn;}
      
      /**
       * function to set position of this atom
       * @param fX x-coordinate of this atom
       * @param fY y-coordinate of this atom
       * @param fZ z-coordinate of this atom
       */
      void setPose(const delphi_real fX, const delphi_real fY, const delphi_real fZ)
      {this->gPose.nX = fX; this->gPose.nY = fY; this->gPose.nZ = fZ;}
      
      /**
       * function to store the rest of line read from .pdb file for this atom
       * @param strAtInfIn string of the line containing the rest info in .pdb for this atom
       */
      void setAtInf(const string strAtInfIn) {this->strAtInf = strAtInfIn;}
      
      /**
       * function to get radius of this atom
       * @return radius of this atom
       */
      delphi_real getRadius() const {return this->fRadius;}
      
      /**
       * function to get charge on this atom
       * @return charge on this atom
       */
      delphi_real getCharge() const {return this->fCharge;}
      
      /**
       * function to get position of this atom
       * @return position of this atom
       */
      SGrid<delphi_real> getPose() const {return this->gPose;}
      
      /**
       * function to get rest info of this atom
       * @return string containing rest info of this atom
       */
      string getAtInf() const {return this->strAtInf;} 
};


#endif // IO_DATATYPE_H_
