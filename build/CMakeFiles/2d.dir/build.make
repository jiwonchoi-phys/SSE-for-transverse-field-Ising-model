# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build"

# Include any dependencies generated for this target.
include CMakeFiles/2d.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/2d.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/2d.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/2d.dir/flags.make

CMakeFiles/2d.dir/2d.cpp.o: CMakeFiles/2d.dir/flags.make
CMakeFiles/2d.dir/2d.cpp.o: /Users/jiwon/OneDrive\ -\ pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/2d.cpp
CMakeFiles/2d.dir/2d.cpp.o: CMakeFiles/2d.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/2d.dir/2d.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/2d.dir/2d.cpp.o -MF CMakeFiles/2d.dir/2d.cpp.o.d -o CMakeFiles/2d.dir/2d.cpp.o -c "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/2d.cpp"

CMakeFiles/2d.dir/2d.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2d.dir/2d.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/2d.cpp" > CMakeFiles/2d.dir/2d.cpp.i

CMakeFiles/2d.dir/2d.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2d.dir/2d.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/2d.cpp" -o CMakeFiles/2d.dir/2d.cpp.s

CMakeFiles/2d.dir/util.cpp.o: CMakeFiles/2d.dir/flags.make
CMakeFiles/2d.dir/util.cpp.o: /Users/jiwon/OneDrive\ -\ pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/util.cpp
CMakeFiles/2d.dir/util.cpp.o: CMakeFiles/2d.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/2d.dir/util.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/2d.dir/util.cpp.o -MF CMakeFiles/2d.dir/util.cpp.o.d -o CMakeFiles/2d.dir/util.cpp.o -c "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/util.cpp"

CMakeFiles/2d.dir/util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2d.dir/util.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/util.cpp" > CMakeFiles/2d.dir/util.cpp.i

CMakeFiles/2d.dir/util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2d.dir/util.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/util.cpp" -o CMakeFiles/2d.dir/util.cpp.s

# Object files for target 2d
2d_OBJECTS = \
"CMakeFiles/2d.dir/2d.cpp.o" \
"CMakeFiles/2d.dir/util.cpp.o"

# External object files for target 2d
2d_EXTERNAL_OBJECTS =

2d: CMakeFiles/2d.dir/2d.cpp.o
2d: CMakeFiles/2d.dir/util.cpp.o
2d: CMakeFiles/2d.dir/build.make
2d: CMakeFiles/2d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable 2d"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/2d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/2d.dir/build: 2d
.PHONY : CMakeFiles/2d.dir/build

CMakeFiles/2d.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/2d.dir/cmake_clean.cmake
.PHONY : CMakeFiles/2d.dir/clean

CMakeFiles/2d.dir/depend:
	cd "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model" "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model" "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build" "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build" "/Users/jiwon/OneDrive - pukyong.ac.kr/Works/SSE-for-transverse-field-Ising-model/build/CMakeFiles/2d.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/2d.dir/depend
