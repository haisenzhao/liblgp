cmake_minimum_required(VERSION 3.5)

# Set the project name
project (pgl)


# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    pgl_functs.hpp
)

SOURCE_GROUP("source" FILES ${SOURCES})

# Add an executable with the above sources
#add_executable(SmartCFS ${SOURCES})

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

set(GlmIncludeDir ${source_dir})


#MESSAGE(STATUS ${CGAL_INCLUDE_DIR} " found")
include_directories(${PROJECT_SOURCE_DIR} ${GlmIncludeDir})

#add_dependencies(${PROJECT_NAME} glm)


