cmake_minimum_required(VERSION 3.21)
project(DCC_Processing_Tool)

set(CMAKE_CXX_STANDARD 17)

add_executable(DCC_Processing_Tool main.cpp src/DCC_Processing/DCC_Processing.h src/DCC_Processing/Processing_Functions.h src/DCC_Kinetic/DCCKinetic.h src/DCC_Kinetic/Kinetic_Functions.h src/DCC_Characterisation/StructureCharacterisation.h src/DCC_Characterisation/TJsLab.h src/DCC_SupportFunctions.h src/DCC_Characterisation/QPsLab.h src/DCC_Characterisation/DislocationsLab.h src/DCC_Writer/DCC_Writer.h src/DCC_Writer/Writer_Functions.h)