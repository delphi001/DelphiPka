/*
 * site_exceptions.h
 *
 *  Created on: Feb 18, 2014
 *      Author: chuan
 */

#ifndef SITE_EXCEPTIONS_H_
#define SITE_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CUnknownGridEngFile : public CException
{
   public:
      CUnknownGridEngFile(const string& strFileName)
      {
         cerr << "data file " << strFileName << " for analytic grid energy not present \n";
      }
};

class CCrgatnError : public CException
{
   public:
      CCrgatnError()
      {
         cerr << "PROBLEM WITH prgiCrgAt (crgatn) \n";
      }
};

class CUnknownPreviousPhiFile : public CException
{
   public:
      CUnknownPreviousPhiFile(const string& strPreviousPhiFile)
      {
         cerr << "THE INPUT PHI FILE " << strPreviousPhiFile << " DOES NOT EXISIT \n";
      }
};

class CUnmatchPotentialMap : public CException
{
   public:
      CUnmatchPotentialMap()
      {
         cerr << "THE TWO POTENTIAL MAPS DO NOT MATCH \n";
      }
};

class CUnknownInFrcFile : public CWarning
{
   public:
      CUnknownInFrcFile(const string & strFrciFile)
      {
         cerr << "THE INPUT FRC FILE " << strFrciFile << " DOES NOT EXISIT ";
         cerr << "(EXITING...)\n";
      }
};

class CNoAtomInfo : public CWarning
{
   public:
      CNoAtomInfo(const string & strFrciFile)
      {
         cerr << "THIS UNFORMATTED FILE " << strFrciFile << " DOES NOT CONTAIN ATOM INFO ";
         cerr << "(ATOM INFO FLAG TRUNED OFF)\n";
      }
};

class CCalcReactForceError : public CWarning
{
   public:
   CCalcReactForceError()
      {
         cerr << "CANNOT CALCULATE REACTION FORCES W/O USING INTERNAL (SELF) COORDINATES ";
         cerr << "(EXITING...)\n";
      }
};

class CSitePhiError : public CWarning
{
   public:
      CSitePhiError()
      {
         cerr << "Something unclear with sitephi array ";
         cerr << "(will be fixed soon...)\n";
      }
};

class CNoIDebMap : public CWarning
{
   public:
      CNoIDebMap()
      {
         cerr << "WRTSIT: THESE SALT CONCENTRATIONS DO NOT HAVE THE BENEFIT OF IDEBMAP ";
         cerr << "(AS YET)\n";
      }
};

class CUntestedPhicon : public CWarning
{
   public:
      CUntestedPhicon()
      {
         cerr << "PHICON: this option has not been tested yet\n";
      }
};

class CNoPotential2CrgConcentrate : public CWarning
{
   public:
      CNoPotential2CrgConcentrate()
      {
         cerr << "CANNOT CONVERT FROM POTENTIALS TO CONCENTRATIONS IF THE IONIC STRENTH IS ZERO\n";
      }
};

class CEmptyPhiMap : public CWarning
{
   public:
      CEmptyPhiMap()
      {
         cerr << "THE REQUESTED OUPUT PHIMAP IS EMPTY.\n";
      }
};

#endif /* SITE_EXCEPTIONS_H_ */
