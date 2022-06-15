cmake_minimum_required(VERSION 3.21)
project(DCC_Processing_Tool)

set(CMAKE_CXX_STANDARD 14)

add_executable(DCC_Processing_Tool MainProcessingTool.cpp src/DCC_SupportFunctions.h)
