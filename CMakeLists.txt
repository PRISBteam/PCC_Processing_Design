cmake_minimum_required(VERSION 3.24)
project(PCC_Processing_Design_tool)

set(CMAKE_CXX_STANDARD 20)

# added -fopenmp
include_directories(/usr/local/Cellar/libomp/15.0.7/include/)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp")

if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}" CACHE STRING "" FORCE)
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "" FORCE)
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
    endif()
endif()

find_package(OpenMP REQUIRED)

add_executable( PCC_Processing_Design_tool
        # Management module
        src/main.cpp
        # Various support functions
        src/DCC_SupportFunctions.h
        # C++  class for PCC and related chains of special elements
        src/DCC_Objects/DCC_Objects.h
        # Processing module
        src/DCC_Processing/DCC_Processing.h
        src/DCC_Processing/Processing_Functions.h
        # Kinetic module
        src/DCC_Kinetic/DCC_Kinetic.h
        src/DCC_Kinetic/Kinetic_Functions.h
        # Characterisation module
        src/DCC_Characterisation/StructureCharacterisation.h
        src/DCC_Characterisation/TJsLab.h
        src/DCC_Characterisation/QPsLab.h
        src/DCC_Characterisation/DislocationsLab.h
        # Writer module
        src/DCC_Writer/DCC_Writer.h
        src/DCC_Writer/Writer_Functions.h
        # Subcomplex module
        src/DCC_Section/Planecut_Functions.h
        # Multiphysics module
        src/DCC_Multiphysics/DCC_Multiphysics.h
        src/DCC_Multiphysics/Multiphysics_Functions.h

        # Simulation Tasks here SHOULD NOT BE !!!
        #src/tasks/task_macrocrack.cpp
        #src/tasks/task_sFacesProcessing.cpp
        )

#target_compile_options(PCC_Processing_Design_tool PRIVATE -Wall ${OpenMP_CXX_FLAGS})
#target_link_libraries(PCC_Processing_Design_tool PRIVATE ${OpenMP_CXX_FLAGS})