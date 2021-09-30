cmake_minimum_required(VERSION 3.5)

# Set the project name
project (pgl)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    pgl_functs.hpp
	main.cpp
)

# Add an executable with the above sources
add_executable(pgl ${SOURCES})

set(GLM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/glm)

include_directories(${PROJECT_SOURCE_DIR} ${GLM_INCLUDE_DIR})


#set_target_properties(pgl PROPERTIES LINKER_LANGUAGE C)