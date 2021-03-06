cmake_minimum_required (VERSION 3.6 FATAL_ERROR)
set(PROJECT_NAME LocalFields)

project (${PROJECT_NAME} LANGUAGES CXX CUDA)
find_package(CUDA QUIET REQUIRED)

set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)

SET(CMAKE_BUILD_TYPE "Release")

set (IGL_DIR "" CACHE PATH "Root directory of libIGL")
set (EIGEN_DIR "" CACHE PATH "Root directory of Eigen")
set (IMGUI_DIR "" CACHE PATH "Library directory of IMGUI")
set (GLFW_INCLUDE_DIR "" CACHE PATH "Include directory of GLFW")
set (GLFW_LIBRARY_DIR "" CACHE PATH "Library directory of GLFW")
set (GLAD_INCLUDE_DIR "" CACHE PATH "Include directory of GLAD")
set (GLAD_LIBRARY_DIR "" CACHE PATH "Library directory of GLAD")

#Defining the system varibles where the include files reside
include_directories(${IGL_DIR}/include/
					${EIGEN_DIR}
					${GLFW_INCLUDE_DIR}/include/
					${GLAD_INCLUDE_DIR}/include/
					${IMGUI_DIR}
)

#Defining the system varibles where the link libraries reside					
link_directories(${GLFW_LIBRARY_DIR}
				 ${GLAD_LIBRARY_DIR}
)	

#Stating what files to include in the project				 
add_executable(${PROJECT_NAME}
	Main_VectorFields.cpp 
	VectorFields.h
	VectorFields.cpp	
	LocalFields.h 
	LocalFields.cpp 
	EigenSolver.h
	EigenSolver.cpp	
	TestSolver.h 
	TestSolver.cpp 
	Utility.h 
	Utility.cpp
	VF_Visualization.cpp
	VF_Testing.cpp
	VF_MatrixConstructions.cpp
)

target_link_libraries(
	${PROJECT_NAME}
	OpenGL32.lib; 			
	glad.lib;
	glfw3.lib;
	cublas.lib; 			cudart.lib; 
	cusolver.lib; 			cusparse.lib;
)	

if(MSVC)	
	add_definitions(-D_USE_MATH_DEFINES)
endif()

#To activate OpenMP
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

set(CMAKE_SUPPRESS_REGENERATION true)