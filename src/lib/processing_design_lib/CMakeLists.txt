project(Processing_Design)

include_directories(../external/Eigen-5.0/src)
include_directories(../external/Spectra-1.2.0)

set(HEADER_FILES
        ../ini/ini_readers.h
        ../ini/ini_materials_reader.h

        PCC_Support_Functions.h
        PCC_Objects.h
        PCC_Measures.h

        PCC_Subcomplex/PCC_Subcomplex.h
        PCC_Subcomplex/functions/subcomplex_cross_section.h
        PCC_Subcomplex/functions/subcomplex_clusters.h

        PCC_Multiphysics/PCC_Multiphysics.cpp
        PCC_Multiphysics/functions/multiphysics_internal_stresses.cpp

        PCC_Processing/PCC_Processing.h
        PCC_Processing/functions/processing_assigned_labelling.h
        PCC_Processing/functions/processing_induced_labelling.h
        PCC_Processing/functions/processing_indexing.h

        PCC_Kinetics/PCC_Kinetics.h
        PCC_Kinetics/functions/corrosion_lab.h
        PCC_Kinetics/functions/irradiation_lab.h

        PCC_Design/PCC_Design.h
        PCC_Design/functions/GeneticAlgorithm.h

        PCC_Writer/PCC_Writer.h
        PCC_Writer/functions/writer_pcc_cells.h
        PCC_Writer/functions/writer_subcomplexes.h
)

set(SOURCE_FILES
        ../ini/ini_readers.cpp
        ../ini/ini_materials_reader.cpp

        PCC_Support_Functions.cpp
        PCC_Objects.cpp
        PCC_Measures.cpp

        PCC_Subcomplex/PCC_Subcomplex.cpp
        PCC_Subcomplex/functions/subcomplex_cross_section.cpp
        PCC_Subcomplex/functions/subcomplex_clusters.cpp

        PCC_Multiphysics/PCC_Multiphysics.h
        PCC_Multiphysics/functions/multiphysics_internal_stresses.h

        PCC_Processing/PCC_Processing.cpp
        PCC_Processing/functions/processing_assigned_labelling.cpp
        PCC_Processing/functions/processing_induced_labelling.cpp
        PCC_Processing/functions/processing_indexing.cpp

        PCC_Kinetics/PCC_Kinetics.cpp
        PCC_Kinetics/functions/corrosion_lab.cpp
        PCC_Kinetics/functions/irradiation_lab.cpp

        PCC_Design/PCC_Design.cpp
        PCC_Design/functions/GeneticAlgorithm.cpp

        PCC_Writer/PCC_Writer.cpp
        PCC_Writer/functions/writer_pcc_cells.cpp
        PCC_Writer/functions/writer_subcomplexes.cpp

        #[===[

    src/lib/processing_design_lib/PCC_Characterisation/PCC_Characterisation.cpp
    src/lib/processing_design_lib/PCC_Characterisation/functions/entropic_analysis.cpp
    src/lib/processing_design_lib/PCC_Characterisation/functions/spectral_analysis.cpp

    src/tasks/cpd_task.cpp
    ]===]#
        )

add_library(processing_design_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})