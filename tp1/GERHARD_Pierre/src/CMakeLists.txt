cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "/usr/share/feel/cmake/modules/")
Find_Package(Feel++)
feelpp_add_application(membrane SRCS membrane.cpp CFG membrane.cfg INCLUDE_IN_ALL)
  
