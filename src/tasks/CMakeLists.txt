project(user_tasks)

set(HEADER_FILES
        user_script.h
        
        )

set(SOURCE_FILES
        user_script.cpp
        
        )

add_library(user_scripts STATIC ${SOURCE_FILES} ${HEADER_FILES})