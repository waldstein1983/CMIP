cmake_minimum_required(VERSION 3.7)
project(CMIP)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
#        ufl_bd.cpp
#        xbcutstk.cxx
        ufl_bc.cpp
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