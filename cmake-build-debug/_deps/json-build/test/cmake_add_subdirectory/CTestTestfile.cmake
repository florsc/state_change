# CMake generated Testfile for 
# Source directory: /home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory
# Build directory: /home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-build/test/cmake_add_subdirectory
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cmake_add_subdirectory_configure "/snap/clion/193/bin/cmake/linux/bin/cmake" "-G" "Ninja" "-Dnlohmann_json_source=/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src" "/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory/project")
set_tests_properties(cmake_add_subdirectory_configure PROPERTIES FIXTURES_SETUP "cmake_add_subdirectory" _BACKTRACE_TRIPLES "/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory/CMakeLists.txt;1;add_test;/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory/CMakeLists.txt;0;")
add_test(cmake_add_subdirectory_build "/snap/clion/193/bin/cmake/linux/bin/cmake" "--build" ".")
set_tests_properties(cmake_add_subdirectory_build PROPERTIES FIXTURES_REQUIRED "cmake_add_subdirectory" _BACKTRACE_TRIPLES "/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory/CMakeLists.txt;7;add_test;/home/florian/CLionProjects/untitled1/cmake-build-debug/_deps/json-src/test/cmake_add_subdirectory/CMakeLists.txt;0;")
