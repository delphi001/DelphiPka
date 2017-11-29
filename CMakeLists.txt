cmake_minimum_required(VERSION 3.8)
project(DelphiPka)

set(CMAKE_VERBOSE_MAKEFILE ON)

if (APPLE)
    set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)
else ()
    set(CMAKE_CXX_COMPILER mpicxx)
endif ()

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-g3")
SET(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(HEADER_FILES
        src/delphi/delphi/delphi_constants.h
        src/delphi/delphi/delphi_data.h
        src/delphi/delphi/delphi_datamarshal.h
        src/delphi/delphi/delphi_exceptions.h
        src/delphi/energy/energy.h
        src/delphi/energy/energy_exceptions.h
        src/delphi/interface/environment.h
        src/delphi/interface/exceptions.h
        src/delphi/interface/interface_exceptions.h
        src/delphi/io/io.h
        src/delphi/io/io_datatype.h
        src/delphi/io/io_exceptions.h
        src/delphi/misc/misc_grid.h
        src/delphi/site/site.h
        src/delphi/site/site_exceptions.h
        src/delphi/solver/solver_exceptions.h
        src/delphi/solver/solver_fastSOR.h
        src/delphi/space/space.h
        src/delphi/space/space_exceptions.h
        src/delphi/space/space_templates.h
        src/delphiPKa/data_store.h
        src/delphiPKa/prime_environment.h
        )

set(SOURCE_FILES
        src/delphi/app/app_primepka.cpp
        src/delphi/delphi/delphi_data_flash.cpp
        src/delphi/delphi/delphi_data_reset.cpp
        src/delphi/delphi/delphi_data_showMap.cpp
        src/delphi/delphi/delphi_data_setMap.cpp
        src/delphi/delphi/delphi_datamarshal_getFunction.cpp
        src/delphi/delphi/delphi_datamarshal_getStatement.cpp
        src/delphi/delphi/delphi_datamarshal_setDefault.cpp
        src/delphi/delphi/delphi_datamarshal_showParameters.cpp
        src/delphi/delphi/delphi_datamarshal_updateParameters.cpp
        src/delphi/interface/exceptions.cpp
        src/delphi/interface/interface_datacontainer.cpp
        src/delphi/interface/interface_datamarshal.cpp
        src/delphi/io/io_epsmap.cpp
        src/delphi/io/io_force.cpp
        src/delphi/io/io_frc.cpp
        src/delphi/io/io_misc.cpp
        src/delphi/io/io_pdb.cpp
        src/delphi/misc/misc_grid_opts.cpp
        src/delphi/misc/misc_interpl.cpp
        src/delphi/misc/misc_timer.cpp
        src/delphi/space/grid_space.cpp
        src/delphi/space/space_crgarr.cpp
        src/delphi/space/space_cube.cpp
        src/delphi/space/space_eps.cpp
        src/delphi/space/space_fft.cpp
        src/delphi/space/space_indver.cpp
        src/delphi/space/space_msrf.cpp
        src/delphi/space/space_run.cpp
        src/delphi/space/space_sas.cpp
        src/delphi/space/space_sclbp.cpp
        src/delphi/space/space_setconvolute.cpp
        src/delphi/space/space_setgaussian.cpp
        src/delphi/space/space_setout.cpp
        src/delphi/space/space_validateInput.cpp
        src/delphi/space/space_vwtms_inc.cpp
        src/delphi/space/space_vwtms.cpp
        src/delphi/solver/solver_bndy_isCoulombBndy.cpp
        src/delphi/solver/solver_bndy_isDipolarBndy.cpp
        src/delphi/solver/solver_bndy_isFocusBndy.cpp
        src/delphi/solver/solver_bndy_setBndy.cpp
        src/delphi/solver/solver_fastSOR_initOddEvenItr.cpp
        src/delphi/solver/solver_fastSOR_itit.cpp
        src/delphi/solver/solver_fastSOR_itrEvenPoints.cpp
        src/delphi/solver/solver_fastSOR_itrOddPoints.cpp
        src/delphi/solver/solver_fastSOR_mkdbsf.cpp
        src/delphi/solver/solver_fastSOR_nitit.cpp
        src/delphi/solver/solver_fastSOR_postItr.cpp
        src/delphi/solver/solver_fastSOR_relfac.cpp
        src/delphi/solver/solver_fastSOR_setcrg.cpp
        src/delphi/solver/solver_fastSOR_validateInput.cpp
        src/delphi/solver/solver_fastSOR_run.cpp
        src/delphi/energy/energy_clb.cpp
        src/delphi/energy/energy_clbmedia.cpp
        src/delphi/energy/energy_clbnonl.cpp
        src/delphi/energy/energy_clbtotal.cpp
        src/delphi/energy/energy_nl.cpp
        src/delphi/energy/energy_react.cpp
        src/delphi/energy/energy_run.cpp
        src/delphi/site/site_expand.cpp
        src/delphi/site/site_phicon.cpp
        src/delphi/site/site_rforce.cpp
        src/delphi/site/site_rforceesp1.cpp
        src/delphi/site/site_tops.cpp
        src/delphi/site/site_writePAnalysis.cpp
        src/delphi/site/site_writePhi.cpp
        src/delphi/site/site_writePhiMap.cpp
        src/delphi/site/site_writePotential_ccp4.cpp
        src/delphi/site/site_writePotential_cube.cpp
        src/delphi/site/site_writePotential_delphi.cpp
        src/delphi/site/site_writePotential_fromPrevious.cpp
        src/delphi/site/site_writePotential_grasp.cpp
        src/delphi/site/site_writePotential_insight.cpp
        src/delphi/site/site_writePotential_zphi.cpp
        src/delphi/site/site_writeSite.cpp
        src/delphiPKa/clustering.cpp
        src/delphiPKa/network.cpp
        src/delphiPKa/energy.cpp
        src/delphiPKa/energy_misc.cpp
        src/delphiPKa/energymap.cpp
        src/delphiPKa/energy_readlog.cpp
        src/delphiPKa/global_type.cpp
        src/delphiPKa/main.cpp
        src/delphiPKa/orbt_type.cpp
        src/delphiPKa/pdb_import.cpp
        src/delphiPKa/place_hydrogen.cpp
        src/delphiPKa/titration2.cpp
        src/delphiPKa/titration_misc.cpp
        src/delphiPKa/readparam.cpp
        src/delphiPKa/topology_import.cpp
        )

add_executable(DelphiPka ${HEADER_FILES} ${SOURCE_FILES})

target_link_libraries(DelphiPka gsl gslcblas)
