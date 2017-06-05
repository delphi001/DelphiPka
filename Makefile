## The Compiler to Use
CC=mpic++

## The Target Executable to Make 
TARGET=delphiPKa

## The Compiler Options to Add
CFLAGS=-O3 -g3 -c -fmessage-length=0 -std=c++11

## The Libraries to Add
LIB=

## The Header Search to Add
INCLUDE=

## The Linker Options to Add
LDFLAGS= -lgsl -lgslcblas

## The Source Code Directories
DELPHIDIR=./src/delphi
DELPHIPKADIR=./src/delphiPKa

## The Source Code to Make
SOURCES=$(DELPHIDIR)/app/app_delphipka.cpp $(DELPHIDIR)/delphi/delphi_data_flash.cpp $(DELPHIDIR)/delphi/delphi_data_reset.cpp $(DELPHIDIR)/delphi/delphi_data_showMap.cpp $(DELPHIDIR)/delphi/delphi_data_setMap.cpp $(DELPHIDIR)/delphi/delphi_datamarshal_getFunction.cpp $(DELPHIDIR)/delphi/delphi_datamarshal_getStatement.cpp $(DELPHIDIR)/delphi/delphi_datamarshal_setDefault.cpp $(DELPHIDIR)/delphi/delphi_datamarshal_showParameters.cpp $(DELPHIDIR)/delphi/delphi_datamarshal_updateParameters.cpp $(DELPHIDIR)/interface/exceptions.cpp $(DELPHIDIR)/interface/interface_datacontainer.cpp $(DELPHIDIR)/interface/interface_datamarshal.cpp $(DELPHIDIR)/io/io_epsmap.cpp $(DELPHIDIR)/io/io_force.cpp $(DELPHIDIR)/io/io_frc.cpp $(DELPHIDIR)/io/io_misc.cpp $(DELPHIDIR)/io/io_pdb.cpp $(DELPHIDIR)/misc/misc_grid_opts.cpp $(DELPHIDIR)/misc/misc_interpl.cpp $(DELPHIDIR)/misc/misc_timer.cpp $(DELPHIDIR)/space/grid_space.cpp $(DELPHIDIR)/space/space_crgarr.cpp $(DELPHIDIR)/space/space_cube.cpp $(DELPHIDIR)/space/space_eps.cpp $(DELPHIDIR)/space/space_indver.cpp $(DELPHIDIR)/space/space_msrf.cpp $(DELPHIDIR)/space/space_run.cpp $(DELPHIDIR)/space/space_sas.cpp $(DELPHIDIR)/space/space_sclbp.cpp $(DELPHIDIR)/space/space_setgaussian.cpp $(DELPHIDIR)/space/space_setout.cpp $(DELPHIDIR)/space/space_validateInput.cpp $(DELPHIDIR)/space/space_vwtms_inc.cpp $(DELPHIDIR)/space/space_vwtms.cpp $(DELPHIDIR)/solver/solver_bndy_isCoulombBndy.cpp $(DELPHIDIR)/solver/solver_bndy_isDipolarBndy.cpp $(DELPHIDIR)/solver/solver_bndy_isFocusBndy.cpp $(DELPHIDIR)/solver/solver_bndy_setBndy.cpp $(DELPHIDIR)/solver/solver_fastSOR_initOddEvenItr.cpp $(DELPHIDIR)/solver/solver_fastSOR_itit.cpp $(DELPHIDIR)/solver/solver_fastSOR_itrEvenPoints.cpp $(DELPHIDIR)/solver/solver_fastSOR_itrOddPoints.cpp $(DELPHIDIR)/solver/solver_fastSOR_mkdbsf.cpp $(DELPHIDIR)/solver/solver_fastSOR_nitit.cpp $(DELPHIDIR)/solver/solver_fastSOR_postItr.cpp $(DELPHIDIR)/solver/solver_fastSOR_relfac.cpp $(DELPHIDIR)/solver/solver_fastSOR_setcrg.cpp $(DELPHIDIR)/solver/solver_fastSOR_validateInput.cpp $(DELPHIDIR)/solver/solver_fastSOR_run.cpp $(DELPHIDIR)/energy/energy_clb.cpp $(DELPHIDIR)/energy/energy_clbmedia.cpp $(DELPHIDIR)/energy/energy_clbnonl.cpp $(DELPHIDIR)/energy/energy_clbtotal.cpp $(DELPHIDIR)/energy/energy_nl.cpp $(DELPHIDIR)/energy/energy_react.cpp $(DELPHIDIR)/energy/energy_run.cpp $(DELPHIDIR)/site/site_expand.cpp $(DELPHIDIR)/site/site_phicon.cpp $(DELPHIDIR)/site/site_rforce.cpp $(DELPHIDIR)/site/site_rforceesp1.cpp $(DELPHIDIR)/site/site_tops.cpp $(DELPHIDIR)/site/site_writePAnalysis.cpp $(DELPHIDIR)/site/site_writePhi.cpp $(DELPHIDIR)/site/site_writePhiMap.cpp $(DELPHIDIR)/site/site_writePotential_ccp4.cpp $(DELPHIDIR)/site/site_writePotential_cube.cpp $(DELPHIDIR)/site/site_writePotential_delphi.cpp $(DELPHIDIR)/site/site_writePotential_fromPrevious.cpp $(DELPHIDIR)/site/site_writePotential_grasp.cpp $(DELPHIDIR)/site/site_writePotential_insight.cpp $(DELPHIDIR)/site/site_writeSite.cpp $(DELPHIPKADIR)/clustering.cpp $(DELPHIPKADIR)/network.cpp $(DELPHIPKADIR)/energy.cpp $(DELPHIPKADIR)/energy_misc.cpp $(DELPHIPKADIR)/energymap.cpp $(DELPHIPKADIR)/energy_readlog.cpp $(DELPHIPKADIR)/global_type.cpp $(DELPHIPKADIR)/main.cpp $(DELPHIPKADIR)/orbt_type.cpp $(DELPHIPKADIR)/pdb_import.cpp $(DELPHIPKADIR)/place_hydrogen.cpp $(DELPHIPKADIR)/titration2.cpp $(DELPHIPKADIR)/titration_misc.cpp $(DELPHIPKADIR)/readparam.cpp $(DELPHIPKADIR)/topology_import.cpp

## The Objects to Generate
OBJECTS=$(SOURCES:.cpp=.o)


## Default rule executed
all: $(TARGET)
	@true    
	
	
## Rule for making the actual target

$(TARGET): $(OBJECTS)
	@echo '	'
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker' 
	$(CC) $(LDFLAGS) $(LIB) $(INCLUDE) $(OBJECTS) -o $@
	@echo 'Finished building target: $@'
	@echo '	'

## Generic compilation rule
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

## Clean Rule	
clean:
	@echo "Clean up .o and delphiPKa"
	@-rm -f $(TARGET) $(OBJECTS)
