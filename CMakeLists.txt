cmake_minimum_required(VERSION 3.22)
project(untitled1)


include(FetchContent)
set(CMAKE_CXX_STANDARD 14)
add_library(simulationLibrary SHARED
        main.cpp)
add_executable(sim main.cpp)


FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json
        GIT_TAG v3.7.3)

FetchContent_GetProperties(json)
if (NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(sim ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(sim nlohmann_json::nlohmann_json)
target_link_libraries(simulationLibrary ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(simulationLibrary nlohmann_json::nlohmann_json)