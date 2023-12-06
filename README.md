This code is header-only Library developed by Haisen Zhao for his research projects.

# Dependency

Depend on [glm](https://github.com/g-truc/glm.git) and [eigen](https://github.com/libigl/eigen.git) but you don't need to install them explicitly.


# Call Liblgp

## Create Project
- Open your Visual Studio and create a new project.![image-create_a_new_project](images/image-create_a_new_project.png)
- Right-click on your project in the "Solution Explorer" and select "Add" "New item..."![add-item](images/add-item.png)

## Call Liblgp in your project
- Download Liblgp\_functs.hpp,RI.hpp,tinyxml2.hpp and local\_libs
- Found Eigen and glm in the local_libs folder and copy Liblgp\_functs.hpp,RI.hpp,tinyxml2.hpp and the Eigen and glm folders under the local folder to the project location![move_file](images/move_file.png)
- Then you can call Liblgp in your code!
- You can run the following code to check if Liblgp can be called

```cpp

#include "iostream"
#include "Liblgp_functs.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

using namespace std;
using namespace Liblgp;

int main(int argc, char* argv[])

{

	std::cerr << "WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	
	system("pause");
	return 0;
}
 
```


# Usage in Cmake

```

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

```
## Cmake option
- When you use Cmake, you can choose third-party library or self-content library.![cmake_option](images/cmake_option.png)

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purposes. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.
