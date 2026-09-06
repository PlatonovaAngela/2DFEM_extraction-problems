cmake_minimum_required(VERSION 3.24)

project(2DFEM LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sdk)

include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR cmake-does-not-exist  
)
FetchContent_MakeAvailable(eigen)

file(GLOB HEADER src/*.h)
file(GLOB SOURCE src/*.cpp)
add_executable(${PROJECT_NAME} ${SOURCE} ${HEADER})
target_include_directories(${PROJECT_NAME} PRIVATE ${eigen_SOURCE_DIR})
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
)