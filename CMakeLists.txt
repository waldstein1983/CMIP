cmake_minimum_required(VERSION 3.7)
project(CMIP)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
#        ufl_bd.cpp
#        xbcutstk.cxx
        ufl_bc_bd.cpp
#        branching.h
#        rounding_heu.h
#        benders_cut.h
#        local_branching.h local_branching.cpp
#        branching.cpp
#        rounding_heu.cpp
#        benders_cut.cpp
#        branch_cut.h branch_cut.cpp branch_cut_lb.h branch_cut_lb.cpp benders_cut_inout.h benders_cut_inout.cpp
        )

set(XPRESS_INC_DIR /opt/xpressmp/include)
set(XPRESS_LINK_DIR /opt/xpressmp/lib)
#
#
include_directories(${XPRESS_INC_DIR})
link_directories(${XPRESS_LINK_DIR})

#add_executable(XpressApplications ${SOURCE_FILES})



add_executable(CMIP ${SOURCE_FILES})

target_link_libraries(CMIP xprb xprl xprnls xprs)