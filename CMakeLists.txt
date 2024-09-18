cmake_minimum_required(VERSION 3.1)
project(GRMSD)
add_library(lib 
        src/main_function.cpp
        src/cube.cpp
        src/local_ICP.cpp
        src/basic_function.cpp
        src/Hungarian.cpp
        src/data_form.cpp)
target_compile_options(lib PUBLIC -O3)
target_compile_features(lib PUBLIC cxx_std_11)


add_executable(GRMSD 
        src/GRMSD.cpp
        )

target_compile_options(GRMSD PUBLIC -O3)
target_compile_features(GRMSD PUBLIC cxx_std_11)

target_link_libraries(GRMSD lib)


add_executable(Substructure 
        src/Substructure.cpp
        )

target_compile_options(Substructure PUBLIC -O3)
target_compile_features(Substructure PUBLIC cxx_std_11)

target_link_libraries(Substructure lib)


add_executable(Substructure_lattice 
        src/Substructure_lattice.cpp
        )

target_compile_options(Substructure_lattice PUBLIC -O3)
target_compile_features(Substructure_lattice PUBLIC cxx_std_11)

target_link_libraries(Substructure_lattice lib)