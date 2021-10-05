cmake_minimum_required(VERSION 3.5)

# Set the project name
project (pgl)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    pgl_functs.hpp
    RI.hpp
    RI.cpp
    tinyxml2.hpp
    tinyxml2.cpp
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
set(GlmIncludeDir ${source_dir})


ExternalProject_Add(
    eigen
    PREFIX ${CMAKE_BINARY_DIR}/eigen
    GIT_REPOSITORY https://github.com/libigl/eigen.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})


target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir} ${EigenIncludeDir})

include_directories(${PROJECT_SOURCE_DIR})

add_dependencies(${PROJECT_NAME} glm)
add_dependencies(${PROJECT_NAME} eigen)

#set_target_properties(pgl PROPERTIES LINKER_LANGUAGE C)