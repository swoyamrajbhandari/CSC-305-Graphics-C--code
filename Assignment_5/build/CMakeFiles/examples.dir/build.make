# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/c/Users/swoya/CSC 305/Assignment_5"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/swoya/CSC 305/Assignment_5/build"

# Include any dependencies generated for this target.
include CMakeFiles/examples.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/examples.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/examples.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/examples.dir/flags.make

CMakeFiles/examples.dir/examples/main.cpp.o: CMakeFiles/examples.dir/flags.make
CMakeFiles/examples.dir/examples/main.cpp.o: ../examples/main.cpp
CMakeFiles/examples.dir/examples/main.cpp.o: CMakeFiles/examples.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/swoya/CSC 305/Assignment_5/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/examples.dir/examples/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/examples.dir/examples/main.cpp.o -MF CMakeFiles/examples.dir/examples/main.cpp.o.d -o CMakeFiles/examples.dir/examples/main.cpp.o -c "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/main.cpp"

CMakeFiles/examples.dir/examples/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/examples.dir/examples/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/main.cpp" > CMakeFiles/examples.dir/examples/main.cpp.i

CMakeFiles/examples.dir/examples/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/examples.dir/examples/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/main.cpp" -o CMakeFiles/examples.dir/examples/main.cpp.s

CMakeFiles/examples.dir/examples/raster.cpp.o: CMakeFiles/examples.dir/flags.make
CMakeFiles/examples.dir/examples/raster.cpp.o: ../examples/raster.cpp
CMakeFiles/examples.dir/examples/raster.cpp.o: CMakeFiles/examples.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/swoya/CSC 305/Assignment_5/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/examples.dir/examples/raster.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/examples.dir/examples/raster.cpp.o -MF CMakeFiles/examples.dir/examples/raster.cpp.o.d -o CMakeFiles/examples.dir/examples/raster.cpp.o -c "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/raster.cpp"

CMakeFiles/examples.dir/examples/raster.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/examples.dir/examples/raster.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/raster.cpp" > CMakeFiles/examples.dir/examples/raster.cpp.i

CMakeFiles/examples.dir/examples/raster.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/examples.dir/examples/raster.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/swoya/CSC 305/Assignment_5/examples/raster.cpp" -o CMakeFiles/examples.dir/examples/raster.cpp.s

# Object files for target examples
examples_OBJECTS = \
"CMakeFiles/examples.dir/examples/main.cpp.o" \
"CMakeFiles/examples.dir/examples/raster.cpp.o"

# External object files for target examples
examples_EXTERNAL_OBJECTS =

examples: CMakeFiles/examples.dir/examples/main.cpp.o
examples: CMakeFiles/examples.dir/examples/raster.cpp.o
examples: CMakeFiles/examples.dir/build.make
examples: CMakeFiles/examples.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/swoya/CSC 305/Assignment_5/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable examples"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/examples.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/examples.dir/build: examples
.PHONY : CMakeFiles/examples.dir/build

CMakeFiles/examples.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/examples.dir/cmake_clean.cmake
.PHONY : CMakeFiles/examples.dir/clean

CMakeFiles/examples.dir/depend:
	cd "/mnt/c/Users/swoya/CSC 305/Assignment_5/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/swoya/CSC 305/Assignment_5" "/mnt/c/Users/swoya/CSC 305/Assignment_5" "/mnt/c/Users/swoya/CSC 305/Assignment_5/build" "/mnt/c/Users/swoya/CSC 305/Assignment_5/build" "/mnt/c/Users/swoya/CSC 305/Assignment_5/build/CMakeFiles/examples.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/examples.dir/depend

