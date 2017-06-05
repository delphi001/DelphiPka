/**
 * @file delphi_datamarshal.h
 * @brief class CDelphiDataMarshal
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the class CDelphiDataMarshal, which inherits from the interface class IDataMarshal.
 * An object of CDelphiDataMarshal defines the variables to be contained in the implemented data container
 * and must be paired with one object of CDelphiData.
 */

#ifndef CDELPHIDATAMARSHAL_H_
#define CDELPHIDATAMARSHAL_H_

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <memory> // std::unique_ptr for CIO

#include "../interface/interface_datacontainer.h"
#include "../interface/interface_datamarshal.h"
#include "../misc/misc_timer.h"
#include "../io/io.h"
#include "delphi_constants.h"
#include "delphi_exceptions.h"

using namespace std;

/**
 * class CDelphiDataMarshal is an implementation of the interface IDataMarshal. It provides not only
 * a particular set of statements and functions allowed in the parameter file, but also the set of variables
 * contained in the private map of class CDelphiDtata and their default values.
 *
 * @note This class must be paired with the class CDelphiData for practical uses. They together define
 * the particular application, namely the delphicpp.
 */
class CDelphiDataMarshal:virtual public IDataMarshal
{
   private:
      //-----------------incl_delphimarshal_getStatement.cpp------------------//
      virtual bool getStatement(string strLineNoSpace);

      /**
       * Function to determine the presence of YES, TRUE, T, ON, FALSE, NO, F, and OFF.
       *
       * @param[in] strArgument The argument where the yes or no are looked for
       * @param[in] strStatement The whole statement. Used to print it to stdout if there is an error.
       * @return    1 if YES/TRUE/T/ON, 0 if FALSE/NO/F/OFF and -1 otherwise
       */
      int yesno(string strArgument, string strStatement);

      //-----------------incl_delphimarshal_getFunction.cpp-------------------//
      virtual bool getFunction(string strLineNoSpace);

      /**
       * Function which removes brackets and non standard characters from the string to obtain the file name
       * or file format
       *
       * @param[in] strArg_UpperCase The argument in upper case.
       * @param[in] strArg_fromInput The argument read from parameter file
       * @return a string in upper case if it indicates the format of the file or a string of the file name
       */
      string getFile_Name_or_Format(string strArg_UpperCase,string strArg_fromInput);

      inline vector<string> getArguments(string strParameters);

      //----------------incl_delphimarshal_showParameters.cpp-----------------//

      /**
       * Function to show read-in parameters
       *
       * part of the outputs in the standard log file
       */
      void showParameters() const;

   public:
      /*
       * delphi valid statements shown in user manual version 5.1
       * see the table of "index of statements and their shorthand" in the manual @ pp.19
       */

      //----------------------- set by Statements ------------------------//

      /**
       * - Long form :
       *    -# AUTOCON
       *    -# AUTOCONVERGENCE
       *    -# AUTOMATICCONVERGENCE
       * - Short form :
       *    -# AUTOC
       * - 2L abbrev. :
       *    -# AC
       * - F95 var.: iautocon
       * - Default : TRUE
       * - Description: \n
       * A flag for automatic convergence. The program by default will automatically calculate the
       * number of iterations needed to attain convergence. It is automatically set if no number of
       * iteration is specified otherwise. \n
       *
       * \note
       * When AUTOC = FLASE, either LINIT or NONIT must be set nonzero so that the solver knows which
       * solver to use and how many iterations to take.
       */
      bool bAutoConverge;

      /**
       * - Long form :
       *    -# BOUNDARYCONDITION
       *    -# BOUNDARYCONDITIONS
       * - Short form :
       *    -# BNDCON
       * - 2L abbrev. :
       *    -# BC
       * - F95 var.: ibctyp
       * - Default : 2(=DIPOLAR)
       * - Description: \n
       * An integer flag specifying the type of boundary condition imposed on the edge of the lattice.
       * Allowed options: \n\n
       * (1) potential is zero. \n\n
       * (2) dipolar. The boundary potentials are approximated by the Debye-Huckel potential of the
       *     equivalent dipole to the molecular charge distribution. Phi is the potential estimated
       *     at a given lattice boundary point, q+ (q-) is the sum of all positive (negative) charges,
       *     and r+(r-) is distance from the point to the center of positive (negative) charge, lambda
       *     is the Debye length. \n\n
       * (3) focusing. The potential map from a previous calculation is read in unit 18, and values
       *     for the potential at the lattice edge are interpolated from this map- clearly the first
       *     map should have been generated with a coarser grid (greater distance between lattice
       *     points) and positioned such that current lattice lies completely within old lattice or
       *     the program will protest. For focusing boundary conditions, the program reads in a
       *     potential map from a previous run, and compares the scale of the focusing map with that
       *     for the current run. If they are the same, it assumes that this is a continuation of a
       *     previous run, and iteration of the potentials contained in the previous potential map is
       *     continued. If the scale is not the same, it checks to ensure that the new lattice lies
       *     completely within the old lattice before interpolating the boundary conditions. \n\n
       * (4) coulombic. They are approximated by the sum of Debye-Huckel potentials of all the
       *     charges. qi is the i'th charge, and ri is the distance from the lattice boundary point
       *     to the charge.
       */
      int iBndyType;

      /**
       * - Long form :
       *    -# BOXFILL
       *    -# PERCENTFILL
       *    -# PERCENTBOXFILL
       * - Short form :
       *    -# PERFIL
       * - 2L abbrev. :
       *    -# PF
       * - F95 var.: perfil
       * - Default : 80.0
       * - Description: \n
       * A percentage of the object longest linear dimension to the lattice linear dimension. This
       * will affect the scale of the lattice (grids/angstrom). The percentage fill of the lattice
       * will depend on the application. A large percentage fill will provide a more detailed mapping
       * of the molecular shape onto the lattice. A perfil less than 20% is not usually necessary or
       * advisable. A very large filling will bring the dielectric boundary of the molecule closer
       * to the lattice edge. This will cause larger errors arising from the boundary potential
       * estimates, which are set to zero or approximated by coulombic/Debye-Huckel-type functions
       * using a uniform solvent dielectric. The error will be minimal for higher salt concentrations
       * or weakly charged molecules. Smaller percentages will increase the accuracy of the boundary
       * conditions, but result in a coarser representation of the molecule. Higher resolution can be
       * achieved more efficiently using focusing. \n
       *
       * \note
       * If the molecule is not centered in the origin of the coordinate system, the perfil reflects
       * the percentage of the system that is actually contained in the lattice. For example, if the
       * maximum dimension of a molecule is 100Å, there is no offset and perfil is 50%, then the box
       * side will be 200Å; but if there is an offset of 20Å in the maximum dimension direction, then
       * the box side will be 280Å. \n\n
       * Scale, grid size and perfil are not independent variables so they cannot all be assigned
       * simultaneously in a single run. In any quantitative calculation, the largest possible scale
       * should be used, preferably greater than 2 grids/angstrom. Without focusing, a perfil of
       * around 50% or 60% is reasonable. For example, if scale is set to 2 and perfil is set to 50,
       * the grid size is calculated automatically given the size of the structure. For larger
       * molecules this could mean a prohibitively large memory requirement. In this case a compromise
       * must be found or focusing could be used. Regardless of grid scale, calculations should be
       * repeated at different scales to assess the size of lattice resolution errors. \n\n
       * A good approach to the calculation could start with a small percentage, say 20%, using
       * Debye-Huckel boundary conditions, and then focus in to say 90% or more, in one (or two)
       * stages, using focusing boundary conditions for the second (and third) runs. It is not
       * necessary for the molecule to lie completely within the grid although then the potential
       * boundary conditions must be generated by focusing. However when calculating solvation
       * energies with box fills of > 100% remember that unexpected results may be obtained since
       * parts of the surface, (and perhaps some charges) are not included in the grid.
       */
      delphi_real fPercentageFill;

      /**
       * - Long form :
       *    -# CHEBIT
       * - Short form :
       *    -# CHEBIT
       * - 2L abbrev. :
       *    -# CI
       * - F95 var.: icheb
       * - Default : FALSE
       * - Description: \n
       * A flag, that if it is true the relaxation parameter for linear convergence process is set
       * equal to 1 (usually not modified from default).
       */
      bool bFixedRelaxParam;

      /**
       * - Long form :
       *    -# CLCSRF
       * - Short form :
       *    -# CLCSRF
       * - 2L abbrev. :
       *    -# CS
       * - F95 var.: isrf
       * - Default : FALSE
       * - Description: \n
       * A flag, that when set to true, outputs a GRASP viewable surface file in the file named
       * grasp.srf.
       */
      bool bOutGraspSurf;

      /**
       * - Long form :
       *    -# CONVERGENCEFRACTION
       * - Short form :
       *    -# CONFRA
       * - 2L abbrev. :
       *    -# CF
       * - F95 var.: icon2
       * - Default : 1
       * - Description: \n
       * A flag that determines the convergence fraction. It decides what fraction of grid points are
       * used in assessing convergence (1=all, 2=half, 5=fifth etc). By default it equals 1 (usually
       * not modified from default).
       */
      int iConvergeFract;

      /**
       * - Long form :
       *    -# CONVERGENCEINTERVAL
       * - Short form :
       *    -# CONINT
       * - 2L abbrev. :
       *    -# CI
       * - F95 var.: icon1
       * - Default : 10
       * - Description: \n
       * A flag that determines at what iteration interval convergence is checked, by default it
       * equals 10.(usually not modified from default) The idea behind this parameter is to allow
       * convergence to be checked less frequently to reduce the amount of time spent.
       */
      int iIterateInterval;

      /**
       * - Long form :
       *    -# EXITUNIFORMDIELECTRIC
       * - Short form :
       *    -# EXITUN
       * - 2L abbrev. :
       *    -# XU
       * - F95 var.: iexun
       * - Default : FALSE
       * - Description: \n
       * A flag to terminate the program if uniform dielectric is present (INDI=EXDI). Usually not
       * modified.
       */
      bool bExitUniformDielect;

      /**
       * - Long form :
       *    -# EXTERIORDIELECTRIC
       *    -# EXTERNALDIELECTRIC
       * - Short form :
       *    -# EXDI
       * - 2L abbrev. :
       *    -# ED
       * - F95 var.: repsout
       * - Default : 80.0
       * - Description: \n
       * The external (solution) dielectric constant. A value of EXDI=1 corresponds to the molecule
       * in vacuum, EXDI=80 to the molecule in water. Depending on the application runs with EXDI
       * equal to either of these values may be used to represent different states in a thermodynamic
       * cycle.
       */
      delphi_real fExDielec;

      /**
       * - Long form :
       *    -# FANCYCHARGE
       *    -# SPHERICALCHARGEDISTRIBUTION
       * - Short form :
       *    -# FCRG
       * - 2L abbrev. :
       *    -# FC
       * - F95 var.: isph
       * - Default : FALSE
       * - Description: \n
       * A flag, normally set to false indicating a linear cubic interpolation of charges to grid
       * points; set to true this turns on a spherical charge interpolation. If an atomic charge does
       * not lie exactly on a grid point, then it must somehow be distributed onto the grid points.
       * If this flag is set false, the standard algorithm is used which distributes a charge to
       * the nearest 8 grid points (quick and simple, see the Proteins paper of Klapper et al.). If
       * this flag is set true, then an algorithm is used which gives a more spherically symmetric
       * charge distribution, although the charge is now spread over a wider region of space. For
       * certain cases this gives higher accuracy for potentials less than 3 grid units from a charge
       * (see Gilson et al. J.Comp. Chem paper), although this point has not been exhaustively explored.
       */
      bool bCrgInterplateType;

      /**
       * - Long form :
       *    -# GRIDCONVERGENCE
       * - Short form :
       *    -# GRDCON
       * - 2L abbrev. :
       *    -# GC
       * - F95 var.: gten
       * - Default : 0.0
       * - Description: \n
       * The value for grid convergence. When set, the criterion used to stop the iterative process is
       * the difference on values of grid energy, this option might slow down the calculation a bit,
       * but provides a very strong criterion.
       */
      delphi_real fGridConverge;

      /**
       * - Long form :
       *    -# GRIDSIZE
       * - Short form :
       *    -# GSIZE
       * - 2L abbrev. :
       *    -# GS
       * - F95 var.: igrid
       * - Default : 0.0
       * - Description: \n
       * An odd integer number of points per side of the cubic lattice, min=5, max=571 (=NGRID,
       * platform dependent). A larger grid size will in general mean a better resolution
       * representation of the molecule on the lattice. This will results in more accurate potentials,
       * but will require more time. The number of iterations required to reach a certain convergence
       * will increase approximately linearly with parameter GS. Since the time per iteration will go
       * up as the cube of this parameter the amount of calculation will thus increase at about the
       * fourth power of GS.
       */
      delphi_integer iGrid;

      /**
       * - Long form :
       *    -# INTERIORDIELECTRIC
       * - Short form :
       *    -# INDI
       * - 2L abbrev. :
       *    -# ID
       * - F95 var.: repsin
       * - Default : 2.0
       * - Description: \n
       * The internal (molecules) dielectric constant. It is used only in single molecule systems for
       * compatibility with the old version. A value of INDI=1 corresponds to a molecule with no
       * polarizability- the state assumed in most molecular mechanics applications. INDI=2 represents
       * a molecule with only electronic polarizability (i.e. assuming no reorientation of fixed
       * dipoles, peptide bonds, etc). A value of 2 is based on the experimentally observed high
       * frequency dielectric behavior of essentially all organic materials. INDI=4-6 represents a
       * process where some small reorganization of molecular dipoles occurs which is not represented
       * explicitly (for example in modeling the effects of site directed mutagenesis experiments,
       * when the structure of the wild type, but not mutant protein is known). According to M.K.
       * Gilson and B. Honig, Biopolymers, 25:2097 (1986) for instance, materials having similar
       * dipole density, dipole moment and flexibility as globular proteins have a dielectric between
       * 4 and 6. In modeling any process where large reorientations of dipoles, or large
       * conformational change occurs, i.e. upon folding or denaturation, using a simple dielectric
       * constant for the molecule would be inappropriate, and the change in conformation should be
       * modeled explicitly.
       */
      delphi_real fInDielec;

      /**
       * - Long form :
       *    -# IONICSTRENGTH/SALT2
       *    -# SALTCONC
       *    -# SALTCONCENTRATION
       * - Short form :
       *    -# SALT/SALT2
       * - 2L abbrev. :
       *    -# IS/S2
       * - F95 var.: conc
       * - Default : 0.0/0.0
       * - Description: \n
       * The concentration of first and second kind of salt,(moles/liter). In the case of a single
       * 1:1 salt, it coincides with ionic strength.
       */
      vector<delphi_real> vctfSalt;

      /**
       * - Long form :
       *    -# IONRADIUS
       * - Short form :
       *    -# IONRAD
       * - 2L abbrev. :
       *    -# IR
       * - F95 var.: exrad
       * - Default : 2.0
       * - Description: \n
       * The thickness of the ion exclusion layer around molecule (Å). IONRAD, in combination with
       * the atomic van der Waals radii in the siz file, determines the regions of space, and hence
       * the lattice points, which are inaccessible to solvent ions. Suggested values is IONRAD = 2.0
       * for sodium chloride. For the purpose of DelPhi, a solvent ion is considered as a point
       * charge, which can approach no closer than its ionic radius, IONRAD, to any atoms van der
       * Waals surface. The ion excluded volume is thus bounded by the contact surface, which is the
       * locus of the ion centres when in van der Waals contact with any accessible atom of the
       * molecule. A zero value for IONRAD will just yield the van der Waals surface. A non zero value
       * of IONRAD will thus introduce a Stern, or ion exclusion layer, around the molecule where the
       * solvent ion concentration will be zero and whose dielectric constant is that of the solvent,
       * EXDI.
       */
      delphi_real fIonRadius;

      /**
       * - Long form :
       *    -# ITERATION
       *    -# ITERATIONS
       *    -# LINEARITERATION
       * - Short form :
       *    -# LINIT
       * - 2L abbrev. :
       *    -# LI
       * - F95 var.: nlit
       * - Default : 0
       * - Description: \n
       * An integer number (> 3) of iterations with linear equation. The convergence behavior of the
       * finite difference procedure is reported in the log file as both the mean and maximum absolute
       * change in potential at the grid points between successive iterations. The latter is probably
       * more important since it puts an upper bound on how much the potential is changing at the grid
       * points. It is suggested that sufficient iterations be performed to give a final maximum
       * change of less than 0.001 kT/e. The number of iterations per se is not important, as long as
       * its sufficient to give the required convergence. The convergence behavior can also be judged
       * from the slope of the semi-log plot of the mean and max changes given in the log file. LINIT
       * is best determined by experience, since the convergence rate depends on several factors.
       * Start with say 100 iterations, and then increase the number of iterations until sufficient.
       * Note that a run can be restarted by using focusing boundary conditions with exactly the same
       * SCALE, PERFIL and ACENTER values (see note 5). Some guidelines are: The number of iterations
       * needed will increase with grid size (GSIZE). It will decrease with decreasing PERFIL, since
       * the potentials converge more rapidly in the solvent. It will decrease with increasing ionic
       * strength. The number is fairly insensitive to the size and number of charges on the molecule.
       */
      int iLinIterateNum;

      /**
       * - Long form :
       *    -# LOGFILECONVERGENCE
       * - Short form :
       *    -# LOGGRP
       * - 2L abbrev. :
       *    -# LG
       * - F95 var.: igraph
       * - Default : FALSE
       * - Description: \n
       * A flag that activates the convergence plot during the run.
       */
      bool bLogGraph;

      /**
       * - Long form :
       *    -# LOGFILEPOTENTIALS
       * - Short form :
       *    -# LOGPOT
       * - 2L abbrev. :
       *    -# LP
       * - F95 var.: ipoten
       * - Default : FALSE
       * - Description: \n
       * A flag that activates the potential listing during the run.
       */
      bool bLogPotential;

      /**
       * - Long form :
       *    -# MAXC
       * - Short form :
       *    -# MAXC
       * - 2L abbrev. :
       *    -# XC
       * - F95 var.: res2
       * - Default : 0.0
       * - Description: \n
       * The convergence threshold value based on maximum change of potential (suggested).
       */
      delphi_real fMaxc;

      /**
       * - Long form :
       *    -# NONLINEARITERATION
       *    -# NONLINEARITERATIONS
       * - Short form :
       *    -# NONIT
       * - 2L abbrev. :
       *    -# NI
       * - F95 var.: nnit
       * - Default : 0
       * - Description: \n
       * An integer number (> = 0) of non-linear iterations. If linear PB equation only is required,
       * NONIT is set to be 0.
       */
      int iNonIterateNum;

      /**
       * - Long form :
       *    -# PERIODICBOUNDARYX/PERIODICBOUNDARYY/PERIODICBOUNDARYZ
       * - Short form :
       *    -# PBX/PBY/PBZ
       * - 2L abbrev. :
       *    -# PX/PY/PZ
       * - F95 var.: iper
       * - Default : FALSE/FALSE/FALSE
       * - Description: \n
       * They are the three logical flags (t/f) for periodic boundary conditions for the x,y,z edges
       * of the lattice respectively. Note that periodic boundary conditions will override other
       * boundary conditions on edges to which they are applied. Periodic boundary conditions can be
       * applied in one or more of the x, y or z directions. When applied, the potential at each
       * periodic lattice boundary point is iterated by supplying its missing neighbor(s) from the
       * corresponding point on the opposite edge of the lattice. This can be used for example to
       * model an infinite length of DNA. Assume that the helical axis of the DNA in the pdb file is
       * aligned along the Z axis. The periodic boundary flags are set to false, false, true, and the
       * percent fill of the box, PERFIL, is adjusted so that an integral number of turns just fill
       * the box in the Z direction. Normal boundary conditions are applied to the X,Y boundaries.
       * By setting two, or three of the boundary flags to true, one can simulate 2 * dimensional
       * or 3 dimensional cubic lattices of molecules.
       *
       * \note
       * iper(1:3) are for periodic boundary conditions on the x,y,z edges and iper(4:6) are for corresponding voltage drop.
       */
      vector<bool> vctbPeriodicBndy;

      /**
       * - Long form :
       *    -# PHICON
       * - Short form :
       *    -# PHICON
       * - 2L abbrev. :
       *    -# N.A.
       * - F95 var.: iconc
       * - Default : FALSE
       * - Description: \n
       * A flag, that maps charge density in a .phi file, with a procedure that is equivalent to the
       * one that saves the potential map. phicon=f produces standard potential output in kT/e
       * (approximately equal to 25.6 mV at 25oC, or to 0.593 kcal/mole of charge). phicon=t will
       * give net solvent ion concentration output in M/l, where for every lattice point inside
       * the molecule the concentration is 0, and the outside concentration is obtained from:
       * (-ionic strength*2*sinh(potential)) or its linearized version if linear PBE is used.
       */
      bool bOutCrgDensity;

      /**
       * - Long form :
       *    -# PROBERADIUS/RADPR2
       * - Short form :
       *    -# PRBRAD/RADPR2
       * - 2L abbrev. :
       *    -# PR/R2
       * - F95 var.: radprb
       * - Default : 1.4/-1.0
       * - Description: \n
       * A radius (Å) of probe molecule that will define solvent accessible surface in the Lee and
       * Richard's sense (relative to the part of the molecule which is internal to an object). In
       * combination with the atomic van der Waals radii in the siz file, PRBRAD determines the
       * regions of space, and hence the lattice points, that are inaccessible to solvent molecules
       * (water). Suggested value is PRBRAD 1.4 for water. To understand how these parameters work,
       * you should be familiar with the concepts of contact and solvent accessible surface, as
       * discussed by Lee and Richards, and by Mike Connolly. For the purpose of DelPhi, any region
       * of space that is accessible to any part of a solvent (water) molecule is considered as
       * having a dielectric of EXDI. A value of zero for PRBRAD used with a siz file containing the
       * standard van der Waals radii values will assign any region of space not inside any atom's
       * van der Waals sphere to the solvent. For more details, please refer to Rocchia et al. J.
       * Comp. Chem. paper.
       */
      vector<delphi_real> vctfProbeRadius;

      /**
       * - Long form :
       *    -# RELAXATIONFACTOR
       * - Short form :
       *    -# RELFAC
       * - 2L abbrev. :
       *    -# RF
       * - F95 var.: uspec
       * - Default : 0.9975
       * - Description: \n
       * The externally assigned value for spectral radius (define spectral radius) (usually not
       * modified from default).
       */
      delphi_real fSpectralRadius;

      /**
       * - Long form :
       *    -# RELPAR
       * - Short form :
       *    -# RELPAR
       * - 2L abbrev. :
       *    -# RR
       * - F95 var.: relpar
       * - Default : 1.0
       * - Description: \n
       * A manually assigned value for relaxation parameter in non-linear iteration convergence process.
       */
      delphi_real fRelaxParam;

      /**
       * - Long form :
       *    -# RMSC
       * - Short form :
       *    -# RMSC
       * - 2L abbrev. :
       *    -# MC
       * - F95 var.: res1
       * - Default : 0.0
       * - Description: \n
       * The convergence threshold value based on maximum change of potential (suggested).
       */
      delphi_real fRmsc;

      /**
       * - Long form :
       *    -# SCALE
       * - Short form :
       *    -# SCALE
       * - 2L abbrev. :
       *    -# SC
       * - F95 var.: scale
       * - Default : 1.2
       * - Description: \n
       * The reciprocal of one grid spacing (grids/angstrom).
       */
      delphi_real fScale;

      /**
       * - Long form :
       *    -# SOLVPB
       * - Short form :
       *    -# SOLVPB
       * - 2L abbrev. :
       *    -# SP
       * - F95 var.: isolv
       * - Default : TRUE
       * - Description: \n
       * A flag, which controls the Poisson-Boltzmann solver. Normally DelPhi will invoke the
       * Poisson-Boltzmann solver but if you are interested in using DelPhi for other things such as
       * calculating surface area or producing a GRASP viewable surface file, you can turn off the
       * solver using this option.
       */
      bool bSolvePB;

      /**
       * - Long form :
       *    -# VAL+1/VAL-1
       * - Short form :
       *    -# VAL+1/VAL-1
       * - 2L abbrev. :
       *    -# +1/-1
       * - F95 var.: ival
       * - Default : 1/1
       * - Description: \n
       * A number > 0, valence of positive (negative) ion constituting salt one.
       */
      vector<int> vctiValence1;

      /**
       * - Long form :
       *    -# VAL+2/VAL-2
       * - Short form :
       *    -# VAL+2/VAL-2
       * - 2L abbrev. :
       *    -# +2/-2
       * - F95 var.: ival2
       * - Default : 0/0
       * - Description: \n
       * A number > 0, valence of positive (negative) ion constituting salt two.
       */
      vector<int> vctiValence2;

      /**
       * - Long form :
       *    -# N.A.
       * - Short form :
       *    -# ATPODS
       * - 2L abbrev. :
       *    -# N.A.
       * - F95 var.: atompotdist
       * - Default : 0.5
       * - Description: \n
       * upper bound of atomic potential for charged atoms (averaged over a spherical surface)
       */
      delphi_real fPotentialUpperBond;

      /**
       * - Long form :
       *    -# N.A.
       * - Short form :
       *    -# TEMPER
       * - 2L abbrev. :
       *    -# N.A.
       * - F95 var.: temperature
       * - Default : 297.3342119
       * - Description: \n
       * temperature in absolute degree
       */
      delphi_real fTemper;

      /**
       * - Long form :
       *    -# N.A.
       * - Short form :
       *    -# VDROPX/VDROPY/VDROPZ
       * - 2L abbrev. :
       *    -# N.A.
       * - F95 var.: vdrop
       * - Default : 0.0/0.0/0.0
       * - Description: \n
       * vdrop%x = “VDROPX”; iper(4)=.true. vdrop%y = “VDROPY”; iper(5)=.true. vdrop%z = “VDROPZ”;
       * iper(6)=.true.
       */
      SGrid<delphi_real> gfPotentialDrop;

      /**
       * - F95 var.: iuspec
       * - Default : FALSE
       * - Description: \n
       * flag for using relaxation factor from prm file
       */
      bool bSpectralRadius;

      /**
       * - F95 var.: imanual
       * - Default : FALSE
       * - Description: \n
       * flag for manual assignment of relaxation parameter
       */
      bool bManualRelaxParam;

      /**
       * - F95 var.: phiintype
       * - Default : 0
       * - Description: \n
       * flag for manual assignment of relaxation parameter
       */
      int iPhiInType;


      //-------------------------- io file names ------------------------//

      /* - IO type : IN
       * - F95 var.: prmnam
       * - Default : fort.10
       * - Description: \n
       * Default extension prm. Contains input parameters.
       */
      //string  strParamFile;      ///< prmnam(fort.10)
      //                           ///< input parameter file

      /**
       * - IO type : IN
       * - F95 var.: siznam
       * - Default : fort.11
       * - Description: \n
       * Default extension siz. List describing the van der Waals radii to be assigned to each
       * atom/residue pdb record type. A sample file is provided together with the code. Note the
       * atom and residue fields ignore case and leading blanks. The residue field may be left blank
       * (wild card), causing a match with the given atom type of any residue. ONLY if the residue
       * field is left blank, the LAST 5 characters of the atom record may be left blank. In this
       * case all atom types beginning with the letter in column 1 will be matched. Records of greater
       * specificity override those of less specificity. Beware of ambiguities like calcium (ca) and
       * alpha carbon! All atoms of an input pdb file must be assigned a radius through the siz
       * file, even if it is 0, or the output will be flagged with a warning.
       */
      string strSizeFile;

      /**
       * - IO type : IN
       * - F95 var.: crgnam
       * - Default : fort.12
       * - Description: \n
       * Default extension crg. List of the atomic charges to be assigned to each
       * atom/residue/number/chain pdb record type. A sample file is provided together with the code.
       * The ascii fields for atom, residue, number and chain ignore case and leading blanks. Any
       * field except the atom name may be left blank and will be treated as a wild card. Records of
       * greater specificity override those of lesser specificity as for the siz file above.
       */
      string strCrgFile;

      /**
       * - IO type : IN
       * - F95 var.: pdbnam
       * - Default : fort.13
       * - Description: \n
       * A Brookhaven protein data bank standard format file containing atom labels and coordinates,
       * or a modified OBJECTFILE. Only records starting with ATOM or HETATM are read; if objects or
       * multi-dielectric option are used, also the keywords MEDIA, OBJECT, CRGDST, DATA are also
       * read. The default extension is pdb. The precise format is essential; using Fortran syntax,
       * (6A1,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3) is used for the atom record. From left
       * to right, the fields contain 'ATOM--' or 'HETATM' atom serial number, atom name, alternate
       * location indicator, residue name, chain identifier, residue sequence number, residue
       * insertion code, x, y, and z coordinates, occupancy, temperature factor, footnote number.
       *
       * \note
       * The program treats the residue number as an ascii string, not as an integer. As a warning
       * to the user, there are many variations, and even outright errors found in the format of pdb
       * files obtained from the web. It would be wise to double-check the contents of a file to
       * save any heartache.
       */
      string strPdbFile;

      /**
       * - IO type : OUT
       * - F95 var.: phinam
       * - Default : fort.14
       * - Description: \n
       * If the flag IBIOS (BIOSYM) is false, then output is in DELPHI format, default extension phi.
       * The output can be either a potential map or a concentration map, with format same as for unit
       * 18 above. The output phi map has the same scale as used in the calculation (i.e, variable)
       * unless format=grasp is specified. The grasp-style phi map format will always interpolate
       * to a 65 x 65 x 65 grid for use in Grasp (or other hardwired display/analysis programs). If
       * the flag IBIOS (BIOSYM) is true, then output is in INSIGHT format, default extension ins.
       * This is an unformatted (binary) file. As it was explained above, the format is provided only
       * for completeness in case that one wants to visualize the file with different than Insight
       * software.
       *
       * \note
       * Note that for grid sizes less than 65, INSIGHT format files will occupy less disk space than
       * the corresponding DELPHI files. ins files are designed as input to a Biosym Corp. stand alone
       * utility called CONTOUR, supplied with INSIGHT Version 2.4. This program will produce contour
       * files for display with INSIGHT.33 \n\n
       * If the flag CUBE is true, then output is in CUBE format (Gaussian Cube). Example:
       * \code
       * Out(phi,file=’phimap.txt’,form=’cube’)
       * \endcode
       * this command creates file ‟phimap.txt‟ in the cube-format. \n\n
       */
      string strPhiFile;

      /**
       * - IO type : IN
       * - F95 var.: frcinam
       * - Default : fort.15
       * - Description: \n
       * Default extension: pdb or frc. List of coordinates where site potentials are output in Unit
       * 16. Format as for Unit 13.
       */
      string strFrciFile;

      /**
       * - IO type : OUT
       * - F95 var.: frcnam
       * - Default : fort.16
       * - Description: \n
       * Default extension frc. A list of potentials and fields at coordinates in pdb file read on
       * unit 15. Format: 12 lines of ascii header information, followed by a variable number of
       * records written as: \n
       * \code
       *       230 format(8G10.3) \n
       *       write(16,230)xo,chrgv,phiv,fx,fy,fz \n
       * \endcode
       * where xo(3) are x ,y ,z coordinates of charge, chrgv is the charge value, phiv is the
       * potential (in kT/e) at that point, and fx, fy, fz are the field components (in kT/e/Å ).
       * The last line of the file is the sum of chrgv*phiv/2 over all the charges in the file. This
       * quantity can be used for calculating solvation and interaction energies.
       */
      string strFrcFile;

      /**
       * - IO type : OUT
       * - F95 var.: epsnam
       * - Default : fort.17
       * - Description: \n
       * Dielectric bit map, default extension: eps. If grid size=65, there are 3*65*65*65 lines
       * joining neighboring grid points, 65*65*65 each in of the x,y,z directions. The midpoint of
       * each line is given a value of 1 if it lies within the solvent accessible volume of the
       * system, 0 if outside. This defines the shape of the molecule and separates the space into
       * different dielectric regions. The format of the output files is described below in case that
       * the user wants to build own software to visualize the map. For compact output purposes the
       * array of INTEGER*4, epsmap(65,65,65,3), is compressed into an INTEGER*2 array, neps(5,65,65),
       * by bit-mapping: the first index of epsmap, range 1-65 is compressed into the first index
       * of neps, range 1-5, where the indices 1-16 go into bits 0-15 of the word with index 1,
       * indices 17-32 -> bits 0-15 of word with index 2 etc. The array neps is then written to an
       * unformatted binary file:
       * \code
       * write (17) imap, scale, oldmid
       * write (17) neps
       * \endcode
       * where imap is an unused integer*4 flag and scale, oldmid(3) are real*4 scaling information
       * as above.
       *
       * \note
       * In the case that the solute is composed of more than one dielectric media, in this release
       * (v.6.1 up to rel. 1.1) the additional information is not included in the fort.17, in order
       * to maintain compatibility with software packages that take it as an input.
       */
      string strEpsFile;

      /**
       * - IO type : IN
       * - F95 var.: phiinam
       * - Default : fort.18
       * - Description: \n
       * Default extension phi, potential map for focusing boundary conditions. Potentials are in
       * kT/e (25.6mV, 0.593 kcal/mole/charge at 25°C). \n\n
       * The format of the file is given below in case that the user wants to adopt the file to its
       * own software. If the users wants to visualize the file with Grasp or Insight, no action
       * should be taken.
       * \code
       *    unformatted (binary file)
       *    character*20 uplbl
       *    character*10 nxtlbl,character*60 toplbl
       *    real*4 phi(65,65,65)
       *    character*16 botlbl
       *    real*4 scale,oldmid(3)
       * \endcode
       * uplbl, nxtlbl, toplbl, botlbl are ascii information. Phi is the 3D array containing values
       * of potential for all the lattice points. Index order is x,y,z. Scale is lattice scale in
       * grid/Å. Oldmid is the x,y,z coordinates in real space (angstroms) of the centre of the
       * lattice: thus the real space coordinates x,y,z of the lattice point for phi(IX,IY,IZ),
       * for the case where IGRID = 65, are:
       * \code
       *   x = (IX - 33)/scale + oldmid(1)
       *   y = (IY - 33)/scale + oldmid(2)
       *   z = (IZ - 33)/scale + oldmid(3)
       * \endcode
       * where 33 = (65+1)/2 is the middle point of the grid.
       */
      string strPhiiFile;

      /**
       * - IO type : OUT
       * - F95 var.: mpdbnam
       * - Default : fort.19
       * - Description: \n
       * If the "modified pdb file" option is activated in a WRITE/OUT function, a logical flag (t/f),
       * iatout, will be set to true and will produce a modified PDB file written on unit 19,
       * containing the: radius and charge assigned to each atom written after the coordinates, in
       * the fields used for occupancy and B factor. It is recommended that this option be set
       * initially so that the user can check that all the radius and charge assignments are correct.
       * An additional check on the charge assignment can be made by looking at the total charge
       * written to the log file.
       */
      string strModifiedPdbFile;

      /**
       * - F95 var.: updbnam
       * - Default : fort.20
       * - Description: \n
       * Output unformatted pdb file
       *
       * \note
       * NOT described in manual!
       */
      string strUnformatPdbFile;

      /**
       * - F95 var.: ufrcnam
       * - Default : fort.21
       * - Description: \n
       * Output unformatted frc file
       *
       * \note
       * NOT described in manual!
       */
      string strUnformatFrcFile;

      /**
       * - F95 var.: srfnam
       * - Default : grasp.srf
       * - Description: \n
       * Output grasp surface file
       *
       * \note
       * NOT described in manual!
       */
      string strGraspFile;

      /**
       * - F95 var.: nrgnam
       * - Default : energy.dat
       * - Description: \n
       * Ouput energy file
       *
       * \note
       * NOT described in manual!
       */
      string strEnergyFile;

      /**
       * - F95 var.: scrgnam
       * - Default : scrg.dat
       * - Description: \n
       * Output surface charge file. Only PDB format is supported.
       *
       * \note
       * NOT described in manual!
       */
      string strScrgFile;

      //----------------------- set by functions ------------------------//
      /**
       * - F95 var.: offset
       * - Set by func.: CENTER or CENT
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       *
       */
      SGrid<delphi_real> gfOffCenter;

      /**
       * - F95 var.: acent
       * - Set by func.: ACENTER or ACENT
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * Used to set oldmid when iacent = .true.
       */
      SGrid<delphi_real> gfAcent;

      /**
       * - F95 var.: iacent
       * - Set by func.: ACENTER or ACENT
       * - Default : FALSE
       * - Description: \n
       * flag for Acent(x,y,z) function
       */
      bool bIsAcent;

      /**
       * - F95 var.: pdbfrm
       * - Set by func.: READ or IN
       * - Default : 10
       * - Description: \n
       * pdb file format, = 0 if unknown format, = 1 if “UN”, = 2 if “MOD”, = 3 if "PQR"
       */
      int iPdbFormatIn;

      /**
       * - F95 var.: ipdbrd
       * - Set by func.: READ or IN
       * - Default : FALSE
       * - Description: \n
       * flag for reading unformatted pdb file
       */
      bool bPdbUnformatIn;

      // set by WRITE or OUT function

      /**
       * - F95 var.: phiwrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing potential map
       */
      bool bPhimapOut;

      /**
       * - F95 var.: phifrm
       * - Set by func.: WRITE or OUT
       * - Default : 0
       * - Description: \n
       * Format of file phinam, = 0 if unknown format,= 1 if “BIOSYS”, = 2 if “GRASP”, = 3 if “CCP4”,
       * = 4 if “DIFF”, = 5 if “CUBE”.
       */
      int iPhiFormatOut;

      /**
       * - F95 var.: ibios
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * Flag for format of output file fort.14. = .false. to output in DELPHI format, = .true. to
       * output in INSIGHT format.
       */
      bool bBiosystemOut;

      /**
       * - F95 var.: ibem
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for saving vertices, normals and triangles in the bem.srf file
       */
      bool bBemSrfOut;

      /**
       * - F95 var.: isite
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for outputting site potential
       */
      bool bSiteOut;

      /**
       * - F95 var.: frcfrm
       * - Set by func.: WRITE or OUT
       * - Default : 0
       * - Description: \n
       * Format of file frcnam, = 0 if unknown format, = 1 if “RC”, = 2 if “R”, = 3 if “UN”
       */
      int iFrcFormatOut;

      /**
       * - F95 var.: epswrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing epsilon map
       */
      bool bEpsOut;

      /**
       * - F95 var.: iatout
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * Flag controlling outputting modified PDB file in UNIT 19
       */
      bool bModPdbOut;

      /**
       * - F95 var.: mpdbfrm
       * - Set by func.: WRITE or OUT
       * - Default : 0
       * - Description: \n
       * Format of file mpdbnam, = 0 if unknown format, = 1 if “PQR”.
       */
      int iModPdbFormatOut;

      /**
       * - F95 var.: ipdbwrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing unformatted pdb file
       */
      bool bUnformatPdbOut;

      /**
       * - F95 var.: ifrcwrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing unformatted frc file
       */
      bool bUnformatFrcOut;

      /**
       * - F95 var.: inrgwrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing energy in unit 42 file
       */
      bool bEngOut;

      /**
       * - F95 var.: iwgcrg
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing grid charge file
       */
      bool bGridCrgOut;

      /**
       * - F95 var.: iacs
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for importing exposed vertices from thehsurf2.dat file
       */
      bool bHsurf2DatOut;

      /**
       * - F95 var.: idbwrt
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for writing DELPHI DB file, containing Debye factor
       */
      bool bDbOut;

      /**
       * - F95 var.: isen
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for calculating surface energy positions and values
       */
      bool bSurfEngOut;

      /**
       * - F95 var.: isch
       * - Set by func.: WRITE or OUT
       * - Default : FALSE
       * - Description: \n
       * flag for outputting surface charge file
       */
      bool bSurfCrgOut;

      /**
       * - F95 var.: scrgfrm
       * - Set by func.: WRITE or OUT
       * - Default : 0
       * - Description: \n
       * Format of file scrgnam, = 0 if unknown format, = 1 if “PDB”.
       */
      int iSurfCrgFormatOut;

      /**
       * - F95 var.: logg
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * flag for calculating grid energy
       */
      bool bGridEng;

      /**
       * - F95 var.: logs
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * flag for calculating solvation energy
       */
      bool bSolvEng;

      /**
       * - F95 var.: logas
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * flag for analytic surface
       */
      bool bAnalySurfEng;

      /**
       * - F95 var.: loga
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * flag for calculating analytic energy
       */
      bool bAnalyEng;

      /**
       * - F95 var.: logions
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * Flag for energy calculation of contribution by the solvent
       */
      bool bIonsEng;

      /**
       * - F95 var.: logc
       * - Set by func.: ENERGY
       * - Default : FALSE
       * - Description: \n
       * flag for calculating coulombic energy
       */
      bool bCoulombEng;

      /**
       * - F95 var.: isita
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report atom info in site function.
       */
      bool bAtomInSite;

      /**
       * - F95 var.: isitq
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report charge in site function
       */
      bool bCrgInSite;

      /**
       * - F95 var.: isitp
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report grid potential info in site function
       */
      bool bGridPotentialInSite;

      /**
       * - F95 var.: isitap
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report atomic potential info in site function
       */
      bool bAtomPotentialInSite;

      /**
       * - F95 var.: isitdeb
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report debye fraction map info in site function
       */
      bool bDebyeFractionInSite;

      /**
       * - F95 var.: isitf
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report field info in site function
       */
      bool bFieldInSite;

      /**
       * - F95 var.: isitr
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report reaction potential info in site function
       */
      bool bReactPotentialInSite;

      /**
       * - F95 var.: isitc
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report coulombic potential info in site function
       */
      bool bCoulombPotentialInSite;

      /**
       * - F95 var.: isitx
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report atomic coordinate information in site function
       */
      bool bAtomCoordInSite;

      /**
       * - F95 var.: isiti
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report salt concentration info in site function
       */
      bool bSaltInSite;

      /**
       * - F95 var.: isitt
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report total potential info in site function
       */
      bool bTotalPotentialInSite;

      /**
       * - F95 var.: isitrf
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report reaction force info in site function
       */
      bool bReactForceInSite;

      /**
       * - F95 var.: isitcf
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report coulombic force info in site function
       */
      bool bCoulombForceInSite;

      /**
       * - F95 var.: isitmd
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag for producing molecular dynamics data for future implementation
       */
      bool bMDInSite;

      /**
       * - F95 var.: isitsf
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report site surface charge and electric field at site's Solvent Accessible Surface
       */
      bool bSurfCrgInSite;

      /**
       * - F95 var.: isittf
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report total force info in site function
       */
      bool bTotalForceInSite;

      /**
       * - F95 var.: isitpot
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to report potential info in site function
       */
      bool bPotentialInSite;

      /**
       * - F95 var.: irea
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag for calculating reaction field and output it in formatted frc file
       */
      bool bReactFieldInFRC;

      /**
       * - F95 var.: iself
       * - Set by func.: SITE
       * - Default : FALSE
       * - Description: \n
       * flag to use processed pdb as frc field for site function activation
       */
      bool bPDB2FRCInSite;

      /**
       * - F95 var.: bufz
       * - Set by func.: BUFFZ
       * - Default : {0,0,0,0,0,0}
       * - Description: \n
       * Defines a box with sides parallel to grid unit vectors that the reaction field energy will
       * then be calculated using ONLY the polarization charges contained in that box. The fixed
       * format is BUFFZ(6i3).
       */
      SExtrema<delphi_integer> eiBuffz;

      /**
       * - F95 var.: ibufz
       * - Set by func.: BUFFZ
       * - Default : FALSE
       * - Description: \n
       * flag activating the BUFFZ feature for reaction field energy
       */
      bool bIsBuffz;

      /**
       * - F95 var.: iTypeSurf
       * - Set by func.: SURFACE
       * - Default : -1
       * - Description: \n
       * flag for surface type
       */
      int iTypeSurf;

      //---------------------------- statements -------------------------//
      array<string,iStatementNum> rgstrStatement_ShortForm, rgstrStatement_2lAbbre;

      //------------------------------ functions ------------------------//
      array<string,iFunctionNum_FullName>   rgstrFunction_FullForm;
      array<string,iFunctionNum_ShortName>  rgstrFunction_ShortForm;

      //------------------------------ DelPhi ---------------------------//
      /**
       * - F95 var.: deblen
       * - Default : 0.0
       * - Description: \n
       * Debye Length value
       */
      delphi_real fDebyeLength;

      /**
       * - F95 var.: epsout
       * - Default : 80.0
       * - Description: \n
       * repsout/epkt
       */
      delphi_real fEpsOut;

      /**
       * - F95 var.: cran
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * range of x-, y- and z- coordinates
       */
      SGrid<delphi_real> gfCoordinateRange;

      /**
       * - F95 var.: pmid
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * system geometric center
       */
      SGrid<delphi_real> gfGeometricCenter;

      /**
       * - F95 var.: oldmid
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * grid box center
       */
      SGrid<delphi_real> gfBoxCenter;

      /**
       * - F95 var.: rionst
       * - Default : 0.0
       * - Description: \n
       * ionic strength
       */
      delphi_real fIonStrength;

      /**
       * - F95 var.: chi1
       * - Default : 0.0
       * - Description: \n
       * Coefficient of x^1 term in Taylor series of the charge concentration.
       */
      delphi_real fTaylorCoeff1;

      /**
       * - F95 var.: chi2
       * - Default : 0.0
       * - Description: \n
       * Coefficient of x^2 term in Taylor series of the charge concentration.
       */
      delphi_real fTaylorCoeff2;

      /**
       * - F95 var.: chi3
       * - Default : 0.0
       * - Description: \n
       * Coefficient of x^3 term in Taylor series of the charge concentration.
       */
      delphi_real fTaylorCoeff3;

      /**
       * - F95 var.: chi4
       * - Default : 0.0
       * - Description: \n
       * Coefficient of x^4 term in Taylor series of the charge concentration.
       */
      delphi_real fTaylorCoeff4;

      /**
       * - F95 var.: chi5
       * - Default : 0.0
       * - Description: \n
       * Coefficient of x^5 term in Taylor series of the charge concentration.
       */
      delphi_real fTaylorCoeff5;

      /**
       * - F95 var.: lognl
       * - Default : FALSE
       * - Description: \n
       * Flag for non linear energy calculation.
       */
      bool bNonlinearEng;

      /**
       * - F95 var.: epkt
       * - Default : 0.0
       * - Description: \n
       * EPKT
       */
      delphi_real fEPKT;

      /**
       * - F95 var.: epsin
       * - Default : 2.0
       * - Description: \n
       * repsin/epkt
       */
      delphi_real fEpsIn;

      /**
       * - F95 var.: ifrcrd
       * - Default : FALSE
       * - Description: \n
       * flag for reading unformatted frc file
       */
      bool bFrcUnformatIn;

      /**
       * - F95 var.: idirectalg
       * - Default : 1
       * - Description: \n
       * Direct mapping of epsilon: (0/1)(n/y)
       */
      int iDirectEpsMap;

      /**
       * - F95 var.: numbmol
       * - Default : 0
       * - Description: \n
       * # of molecules
       */
      delphi_integer iMoleculeNum;

      /**
       * - F95 var.: rdmx
       * - Default : 0.01
       * - Description: \n
       * max among the radii
       */
      delphi_real fMaxRadius;

      /**
       * - F95 var.: uniformdiel
       * - Default : TRUE
       * - Description: \n
       * true if the dielectric in the system is uniform
       */
      bool bUniformDielec;

      /**
       * - F95 var.: limobject(Nobject)
       * - Default : AUTOMATIC
       * - Description: \n
       * contains extreme values of each object, for a molecule it has extreme but without radii
       */
      vector< SExtrema<delphi_real> > vctefExtrema;

      /**
       * - F95 var.: xn1(natom)
       * - Default : AUTOMATIC
       * - Description: \n
       * atom coordinates in angstroms
       */
      vector< SGrid<delphi_real> > vctgfAtomCoordA;

      /**
       * - F95 var.: xn2(natom)
       * - Default : AUTOMATIC
       * - Description: \n
       * atom coordinates in grid units
       */
      vector< SGrid<delphi_real> > vctgfAtomCoordG;

      //-------------------------------- IO -----------------------------//

      /**
       * - F95 var.: nmedia
       * - Default : 1
       * - Description: \n
       * # of media
       */
      delphi_integer iMediaNum;

      /**
       * - F95 var.: nobject
       * - Default : 1
       * - Description: \n
       * # of objects
       */
      delphi_integer iObjectNum;

      /**
       * - F95 var.: natom
       * - Default : 0
       * - Description: \n
       * # of atoms
       */
      delphi_integer iAtomNum;

      /**
       * - F95 var.: resnummax
       * - Default : 0
       * - Description: \n
       * maximum residue number
       */
      delphi_integer iResidueNum;

      /**
       * - F95 var.: ionlymol
       * - Default : TRUE
       * - Description: \n
       * true if there are only molecules in the system (no objects)
       */
      bool bOnlyMolecule;

      /**
       * - F95 var.: delphipdb(natom)
       * - Default : AUTOMATIC
       * - Description: \n
       * array of structure to store info read from pdb file
       */
      vector<CAtomPdb> vctapAtomPdb;

      /**
       * - F95 var.: medeps(0:nmedia)
       * - Default : AUTOMATIC
       * - Description: \n
       * vector containing correspondence media<->epsilon/epkt
       */
      vector<delphi_real> vctfMediaEps;

      /**
       * - F95 var.: dataobject(nobject,2)
       * - Default : AUTOMATIC
       * - Description: \n
       * vector containing string with object data, and pre-elab data changed it to vctstrObject(2*nobjectmax)
       */
      vector<string> vctstrObject;

      /**
       * - F95 var.: iatmmed(Natom+Nobjectmax)
       * - Default : AUTOMATIC
       * - Description: \n
       * vector containing internal media-number per atom and object
       */
      vector<delphi_integer> vctiAtomMediaNum;

      //------------------------------ Surface --------------------------//

      /**
       * - F95 var.: nqass
       * - Default : 0
       * - Description: \n
       * number of assigned charges
       */
      delphi_integer iCrgGridNum;

      /**
       * - F95 var.: qnet
       * - Default : 0.0
       * - Description: \n
       * net assigned charge
       */
      delphi_real fNetCrg;

      /**
       * - F95 var.: qmin
       * - Default : 0.0
       * - Description: \n
       * assigned negative charge
       */
      delphi_real fMinusCrg;

      /**
       * - F95 var.: qplus
       * - Default : 0.0
       * - Description: \n
       * assigned positive charge
       */
      delphi_real fPlusCrg;

      /**
       * - F95 var.: cqplus
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * Center of assigned positive charge
       */
      SGrid<delphi_real> gfPlusCrgCenter;

      /**
       * - F95 var.: cqmin
       * - Default : {0.0,0.0,0.0}
       * - Description: \n
       * Center of assigned negative charge
       */
      SGrid<delphi_real> gfMinusCrgCenter;

      /**
       * - F95 var.: cmin
       * - Default : {6000.0,6000.0,6000.0}
       * - Description: \n
       * minimal x-, y- and z- coordinates
       */
      SGrid<delphi_real> gfMinCoordinate;

      /**
       * - F95 var.: cmax
       * - Default : {-6000.0,-6000.0,-6000.0}
       * - Description: \n
       * maximal x-, y- and z- coordinates
       */
      SGrid<delphi_real> gfMaxCoordinate;

      /**
       * - F95 var.: ibnum
       * - Default : 0
       * - Description: \n
       * # of boundary grid points
       */
      delphi_integer iBndyGridNum;

      /**
       * - F95 var.: iepsmp(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * a listing of boundary elements 3D eps map, used to constrcut db array. (Can't get rid of for
       * calculating the nonlinear energy)
       */
      vector< SGrid<delphi_integer> > vctgiEpsMap;

      /**
       * - F95 var.: idebmap(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * logical 3D array for assigning dielectric constants used for the molecular surface scaling
       */
      vector<bool> vctbDielecMap;

      /**
       * - F95 var.: ibgrd(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       * boundary grids
       */
      vector< SGrid<delphi_integer> > vctgiBndyGrid;

      /**
       * - F95 var.: nqgrd
       * - Default : 0
       * - Description: \n
       * number of charges which will be charging the grid
       */
      delphi_integer iCrg2GridNum;

      /**
       * - F95 var.: chrgv2(nqgrd)
       * - Default : AUTOMATIC
       * - Description: \n
       * charges which will be charging the grid
       */
      vector< SGridValue<delphi_real> > vctgvfCrg2Grid;

      /**
       * - F95 var.: nqgrdtonqass(nqgrd)
       * - Default : AUTOMATIC
       * - Description: \n
       * nqgrdtonqass maps ic2-th charge internal atmeps6 to ic1-th general charge
       */
      vector<delphi_integer> vctiCrg2GridMap;

      /**
       * - F95 var.: atmcrg(nqass)
       * - Default : AUTOMATIC
       * - Description: \n
       * atmcrg contains grid positions of all charges AND the charge in the 4th field
       */
      vector< SGridValue<delphi_real> > vctgvfAtomCrg;

      /**
       * - F95 var.: chgpos(nqass)
       * - Default : AUTOMATIC
       * - Description: \n
       * charge position in angstroms
       */
      vector< SGrid<delphi_real> > vctgfCrgPoseA;

      /**
       * - F95 var.: scspos(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       * position in angstroms of induced surface charges
       */
      vector< SGrid<delphi_real> > vctgfSurfCrgA;

      /**
       * - F95 var.: crgatn(nqass)
       * - Default : AUTOMATIC
       * - Description: \n
       *
       */
      vector<delphi_integer> vctiCrgAt;

      /**
       * - F95 var.: atsurf(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       *
       */
      vector<delphi_integer> vctiAtSurf;

      /**
       * - F95 var.: atndx(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       *
       */
      vector<delphi_integer> vctiAtNdx;

      /**
       * - F95 var.: scsnor(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       *
       */
      vector< SGrid<delphi_real> > vctgfSurfCrgE;

      /**
       * - F95 var.: atmeps(nqass)
       * - Default : AUTOMATIC
       * - Description: \n
       *
       */
      vector<delphi_real> vctfAtomEps;

      //------------------------ Gaussian & MEMPOT -----------------------//
      //Lin Li: this is for Gaussian and MEMPOT options
       /**
       * - F95 var.: cutoff
       * - Default : 0
       * - Description: \n
       *
       */
      float fCutoff;

       /**
       * - F95 var.: sigma
       * - Default : 0
       * - Description: \n
       *
       */
      float fSigma;

      /**
       * - F95 var.: inhomo
       * - Default : 0
       * - Description: \n
       *
       */
      int iInhomo;

       /**
       * - F95 var.: srfcut
       * - Default : 0
       * - Description: \n
       *
       */
      float fSrfcut;

      /**
       * - F95 var.: Gaussian
       * - Default : 0
       * - Description: \n
       *
       */
      int iGaussian;

      /**
       * - F95 var.: gepsmp(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * a listing of boundary elements 3D eps map, used to constrcut db array. (Can't get rid of for
       * calculating the nonlinear energy)
       */
      vector< SGrid<delphi_real> > vctgfGepsMap;

      /**
       * - F95 var.: gepsmp2(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * a listing of boundary elements 3D eps map, used to constrcut db array. (Can't get rid of for
       * calculating the nonlinear energy)
       */
      vector< SGrid<delphi_real> > vctgfGepsMap2;

       /**
       * - F95 var.: ergsgaussian
       * - Default : 0
       * - Description: \n
       *
       */
      delphi_real fErgsgaussian;


       /**
       * - F95 var.: radipz
       * - Default : 0
       * - Description: \n
       *
       */
      float fRadipz;

      //------------------------------ Solver ---------------------------//

      /**
       * - F95 var.: icount2b
       * - Default : 0
       * - Description: \n
       * used for realigning idpos and db,compressing to contiguous space
       */
      delphi_integer iDielecBndySum;

      /**
       * - F95 var.: icount1b
       * - Default : 0
       * - Description: \n
       * total charged grid points
       */
      delphi_integer iCrgedGridSum;

      /**
       * - F95 var.: gchrg(icount1b)
       * - Default : AUTOMATIC
       * - Description: \n
       * fractional charge in electron units assigned to each grid point
       */
      vector<delphi_real> vctfGridCrg;

      /**
       * - F95 var.: gchrgp(icount1b)
       * - Default : AUTOMATIC
       * - Description: \n
       * position of each such charge on the grid
       */
      vector< SGrid<delphi_integer> > vctgiGridCrgPose;

      /**
       * - F95 var.: ibc
       * - Default : 0
       * - Description: \n
       * number of charged boundary grid points
       */
      delphi_integer iCrgBdyGrid;

      /**
       * - F95 var.: cgbp(ibc)
       * - Default : AUTOMATIC
       * - Description: \n
       * information on the charged boundary grid points
       */
      vector<SDoubleGridValue> vctdgvCrgBndyGrid;

      /**
       * - F95 var.: phimap(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * 3D potential map
       */
      vector<delphi_real> vctfPhiMap;

      /**
       * - F95 var.: phimap_pre(igrid,igrid,igrid)
       * - Default : AUTOMATIC
       * - Description: \n
       * 3D potential map
       */
      vector<delphi_real> vctfPhiMap_Pre;

      //------------------------------ Energy ---------------------------//

      /**
       * - F95 var.: schrg(ibnum)
       * - Default : AUTOMATIC
       * - Description: \n
       * the induced surface charges in electrons
       */
      vector<delphi_real> vctfSurfCrgE;

      /**
       * - F95 var.: test_ergg
       * - Default : 0.0
       * - Description: \n
       * total grid energy
       */
      delphi_real fEngGrid;

      /**
       * - F95 var.: test_ergc
       * - Default : 0.0
       * - Description: \n
       * coulombic energy
       */
      delphi_real fEngCoul;

      /**
       * - F95 var.: test_ergs
       * - Default : 0.0
       * - Description: \n
       * corrected reaction field energy
       */
      delphi_real fEngCorrect;

      /**
       * - F95 var.: test_ergr
       * - Default : 0.0
       * - Description: \n
       * total reaction field energy
       */
      delphi_real fEngReact;

      /**
       * - F95 var.: test_ergions
       * - Default : 0.0
       * - Description: \n
       * total ionic direct contribution
       */
      delphi_real fEngIons;

      //------------------------- NOT TO BE MAPPED ----------------------//

      /**
       * - F95 var.: centnam
       * - Default : fort.27
       * - Description: \n
       * Site coordinates file. List of coordinates where site potentials are output in Unit 16.
       * Format as for Unit 13.
       */
      string strCentFile;

      /**
       * - F95 var.: rmaxdim
       * - Default : 0.0
       * - Description: \n
       * largest dimension
       */
      delphi_real   fMaxDimension;

      shared_ptr<CTimer> pTimer;   //pointer to an object of CTimer class

      /**
       * - Default : empty
       * - Description: \n
       * String vector of PDB or PQR information communicated between DelPhi interface and PrimePKA in I/O module
       */
      vector<string> strCommPDB;

#ifdef PRIME
      /**
       * - Default : empty
       * - Description: \n
       * String vector of FRC Input communicated between DelPhi interface and PrimePKA in Site function
       */
      vector<string> strCommFRCIn;
#endif
      /**
       * - Default : false
       * - Description: \n
       * The switch to control SiteWrite() for FRC Input communicated between DelPhi interface and PrimePKA
       */
      bool bCommFRCIn;

      /**
       * constructor to generate regular stand-alone executable delphicpp.
       *
       * @param[in] argc   Number of parameters in the command line
       * @param[in] argv[] paraemters in the command line
       * @param[in] pt     pointer to an object of class CTimer to report execution time
       */
      CDelphiDataMarshal(int argc,char* argv[],shared_ptr<CTimer> pt):IDataMarshal(argc,argv),pTimer(pt)
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*             CDelphiDataMarshal is constructed                *\n";
         cout << "****************************************************************\n";
#endif

         setDefault();
      };

      /**
       * constructor to allow delphicpp to be compiled with mcce in order to avoid intensive IO operations.
       *
       * @param[in] mcce_data A pointer to the interface struct SMCCE
       * @param[in] pt        pointer to an object of class CTimer to report execution time
       */
      CDelphiDataMarshal(shared_ptr<CTimer> pt):IDataMarshal(),pTimer(pt)
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*             CDelphiDataMarshal is constructed                *\n";
         cout << "****************************************************************\n";
#endif

         setDefault();
      };

      /**
       * destructor
       */
      ~CDelphiDataMarshal()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                CDelphiDataMarshal is destroyed               *\n";
         cout << "****************************************************************\n";
#endif
      };

      /**
       * set default values for all variables contained in data container
       */
      void setDefault();

      /**
       * function implementing post-reading updates of parameters
       */
      virtual void updateParameters();
};

#endif // CDELPHIDATAMARSHAL_H_
