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



include(ExternalProject)
ExternalProject_Add(
    glm
    PREFIX ${CMAKE_BINARY_DIR}/glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(glm source_dir)

set(GlmIncludeDir ${source_dir}/glm)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir})


include_directories(${PROJECT_SOURCE_DIR})

add_dependencies(${PROJECT_NAME} glm)

#set_target_properties(pgl PROPERTIES LINKER_LANGUAGE C)