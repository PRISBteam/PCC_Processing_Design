cmake_minimum_required(VERSION 3.22)
project(PCC_Processing_Design)

set(CMAKE_CXX_STANDARD 17)

add_executable( PCC_Processing_Design
        # 1 # Simulation environment and the management module
        src/main.cpp
        src/lib/ini/ini_readers.cpp
        src/lib/PCC_Support_Functions.cpp
        src/lib/PCC_Objects.cpp
        src/lib/PCC_Measures.cpp
        src/lib/PCC_Processing/PCC_Processing.cpp
        src/lib/PCC_Processing/functions/Processing_Assignment_functions.cpp
        src/lib/PCC_Writer/PCC_Writer.cpp
        src/lib/PCC_Writer/functions/Writer_functions.cpp
        # 2 # Various support functions
        # 3 # C++ classes and related sets of special elements
        # 4 # Processing module
        # 9 # Writer module
        # 10 # Simulation Tasks
        )