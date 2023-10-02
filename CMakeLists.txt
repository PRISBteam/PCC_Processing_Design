cmake_minimum_required(VERSION 3.22)
project(PCC_Processing_Design)

set(CMAKE_CXX_STANDARD 17)

add_executable( PCC_Processing_Design
        # 1 Simulation environment and the management module
        src/main.cpp
        # 2 Various support functions  # 3 C++ classes and related sets of special elements  # 4 Processing module
        # 5 Kinetic module # 6 Characterisation module # 7 Subcomplex module  # 8 Multiphysics module
        # 9 Writer module # Simulation Tasks
        )