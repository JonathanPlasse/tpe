project(tpe)

cmake_minimum_required(VERSION 2.8)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    SET(GCC_COMPILE_FLAGS "-DLINUX")
ELSE () 
   SET(GCC_COMPILE_FLAGS "") 
   
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

add_definitions(${GCC_COMPILE_FLAGS})
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}" )

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

