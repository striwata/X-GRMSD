# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.21.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.21.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build

# Include any dependencies generated for this target.
include CMakeFiles/GRMSD.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/GRMSD.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/GRMSD.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GRMSD.dir/flags.make

CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o: CMakeFiles/GRMSD.dir/flags.make
CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o: ../src/GRMSD.cpp
CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o: CMakeFiles/GRMSD.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o -MF CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o.d -o CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o -c /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/src/GRMSD.cpp

CMakeFiles/GRMSD.dir/src/GRMSD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GRMSD.dir/src/GRMSD.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/src/GRMSD.cpp > CMakeFiles/GRMSD.dir/src/GRMSD.cpp.i

CMakeFiles/GRMSD.dir/src/GRMSD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GRMSD.dir/src/GRMSD.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/src/GRMSD.cpp -o CMakeFiles/GRMSD.dir/src/GRMSD.cpp.s

# Object files for target GRMSD
GRMSD_OBJECTS = \
"CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o"

# External object files for target GRMSD
GRMSD_EXTERNAL_OBJECTS =

GRMSD: CMakeFiles/GRMSD.dir/src/GRMSD.cpp.o
GRMSD: CMakeFiles/GRMSD.dir/build.make
GRMSD: liblib.a
GRMSD: CMakeFiles/GRMSD.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable GRMSD"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GRMSD.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GRMSD.dir/build: GRMSD
.PHONY : CMakeFiles/GRMSD.dir/build

CMakeFiles/GRMSD.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GRMSD.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GRMSD.dir/clean

CMakeFiles/GRMSD.dir/depend:
	cd /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build /Users/nabikatomonohiro/thesis/c++_GoICP/G-RMSD/build/CMakeFiles/GRMSD.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GRMSD.dir/depend
