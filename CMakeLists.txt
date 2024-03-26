cmake_minimum_required(VERSION 3.20)
project(PCC_Processing_Design)

set(CMAKE_CXX_STANDARD 17)

add_executable( PCC_Processing_Design
        src/main.cpp
        src/lib/ini/ini_readers.cpp

        src/lib/PCC_Support_Functions.cpp
        src/lib/PCC_Objects.cpp
        src/lib/PCC_Measures.cpp

        #        src/lib/PCC_Section/PCC_Subcomplex.cpp
        #       src/lib/PCC_Section/functions/Subcomplex_Planecut_functions.cpp
        src/lib/PCC_Processing/PCC_Processing.cpp
        src/lib/PCC_Processing/functions/Processing_Assignment_functions.cpp

        src/lib/PCC_Writer/PCC_Writer.cpp
        src/lib/PCC_Writer/functions/Writer_functions.cpp
        )