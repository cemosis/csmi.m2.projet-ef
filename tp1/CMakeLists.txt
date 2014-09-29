cmake_minimum_required(VERSION 2.8) 
project(lap2)
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} ) 
  FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH
    $ENV{FEELPP_DIR}/share/feel/cmake/modules /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
  if ( FEELPP_CMAKE_MODULES ) 
    set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
  else() 
    message(FATAL_ERROR      "Feel++ does not seem to have\ been installed on this platform")
  endif()
  Find_Package(Feel++) 
endif()
feelpp_add_application(lap2 SRCS lap2.cpp)
  
