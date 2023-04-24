cmake_minimum_required(VERSION 3.24)
project(PCC_Processing_Design)

set(CMAKE_CXX_STANDARD 20)

# added -fopenmp
include_directories(/usr/local/Cellar/libomp/15.0.7/include/)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp")

if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        #message( "message to display")
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

add_executable( PCC_Processing_Design
        # 1 Simulation environment and management module
        src/main.cpp
        # 2 Various support functions
#        src/lib/PCC_SupportFunctions.h
        # 3 C++  class for PCC and related sets of special elements
#        src/lib/PCC_Objects.h
        # 4 Processing module
#        src/PCC_Processing/PCC_Processing.h
#        src/PCC_Processing/Processing_Functions.h
        # 5 Kinetic module
#        src/PCC_Kinetic/PCC_Kinetic.h
#        src/PCC_Kinetic/Kinetic_Functions.h
        # 6 Characterisation module
#        src/PCC_Characterisation/StructureCharacterisation.h
#        src/PCC_Characterisation/TJsLab.h
#        src/PCC_Characterisation/QPsLab.h
#        src/PCC_Characterisation/DislocationsLab.h
        # 7 Subcomplex module
#        src/PCC_Section/Planecut_Functions.h
        # 8 Multiphysics module
#        src/PCC_Multiphysics/PCC_Multiphysics.h
#        src/PCC_Multiphysics/Multiphysics_Functions.h
        # 9 Writer module
#        src/PCC_Writer/PCC_Writer.h
#        src/PCC_Writer/Writer_Functions.h

        # Simulation Tasks here SHOULD NOT BE !!!
        #src/tasks/task_macrocrack.cpp
        #src/tasks/task_sFacesProcessing.cpp
        #src/tasks/task_StripsProcessingDesign.cpp

        #src/lib/PCC_Processing/functions/Procesing_Induced_functions.h src/lib/PCC_Characterisation/Characterisation_Indexing_functions.h src/lib/PCC_Processing/functions/Processing_Imposed_functions.h
        src/lib/measures.h)

#target_compile_options(PCC_Processing_Design PRIVATE -Wall ${OpenMP_CXX_FLAGS})
#target_link_libraries(PCC_Processing_Design PRIVATE ${OpenMP_CXX_FLAGS})