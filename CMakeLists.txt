cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)

# Set the project name
project (liblgp)

#add _CRT_SECURE_NO_WARNINGS
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

#define option to use local libs
option(USE_LOCAL_LIBS "Use local library" OFF)


# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    Liblgp_functs.hpp
    RI.hpp
    tinyxml2.hpp
    main.cpp
)
# Add an executable with the above sources
add_executable(liblgp ${SOURCES})

if(USE_LOCAL_LIBS)
   set(GlmIncludeDir ${PROJECT_SOURCE_DIR}/local_libs)
   set(EigenIncludeDir ${PROJECT_SOURCE_DIR}/local_libs)
else()
   include(ExternalProject)
   ExternalProject_Add(
       glm
       PREFIX ${CMAKE_BINARY_DIR}/third_party/glm
       GIT_REPOSITORY https://github.com/g-truc/glm.git
       CONFIGURE_COMMAND ""
   	   UPDATE_DISCONNECTED 1
       BUILD_COMMAND ""
       INSTALL_COMMAND ""
       LOG_DOWNLOAD ON
       )
   ExternalProject_Get_Property(glm source_dir)
   set(GlmIncludeDir ${source_dir})


   ExternalProject_Add(
       eigen
       PREFIX ${CMAKE_BINARY_DIR}/third_party/eigen
       GIT_REPOSITORY https://github.com/libigl/eigen.git
       CONFIGURE_COMMAND ""
	   UPDATE_DISCONNECTED 1
       BUILD_COMMAND ""
       INSTALL_COMMAND ""
       LOG_DOWNLOAD ON
       )
   ExternalProject_Get_Property(eigen source_dir)
   set(EigenIncludeDir ${source_dir})
endif()



target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR} PRIVATE ${GlmIncludeDir} ${EigenIncludeDir})

include_directories(${PROJECT_SOURCE_DIR})

if(NOT USE_LOCAL_LIBS)
   add_dependencies(${PROJECT_NAME} glm)
   add_dependencies(${PROJECT_NAME} eigen)
endif()

#set_target_properties(liblgp PROPERTIES LINKER_LANGUAGE C)
