# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/ylang/DATA/I-mode/BOUT-dev

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/ylang/DATA/I-mode/BOUT-dev

# Include any dependencies generated for this target.
include externalpackages/PVODE/CMakeFiles/pvpre.dir/depend.make

# Include the progress variables for this target.
include externalpackages/PVODE/CMakeFiles/pvpre.dir/progress.make

# Include the compile flags for this target's objects.
include externalpackages/PVODE/CMakeFiles/pvpre.dir/flags.make

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o: externalpackages/PVODE/CMakeFiles/pvpre.dir/flags.make
externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o: externalpackages/PVODE/precon/pvbbdpre.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ylang/DATA/I-mode/BOUT-dev/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o -c /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/pvbbdpre.cpp

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.i"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/pvbbdpre.cpp > CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.i

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.s"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/pvbbdpre.cpp -o CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.s

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.o: externalpackages/PVODE/CMakeFiles/pvpre.dir/flags.make
externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.o: externalpackages/PVODE/precon/band.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ylang/DATA/I-mode/BOUT-dev/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.o"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pvpre.dir/precon/band.cpp.o -c /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/band.cpp

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pvpre.dir/precon/band.cpp.i"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/band.cpp > CMakeFiles/pvpre.dir/precon/band.cpp.i

externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pvpre.dir/precon/band.cpp.s"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/precon/band.cpp -o CMakeFiles/pvpre.dir/precon/band.cpp.s

# Object files for target pvpre
pvpre_OBJECTS = \
"CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o" \
"CMakeFiles/pvpre.dir/precon/band.cpp.o"

# External object files for target pvpre
pvpre_EXTERNAL_OBJECTS =

externalpackages/PVODE/libpvpre.a: externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/pvbbdpre.cpp.o
externalpackages/PVODE/libpvpre.a: externalpackages/PVODE/CMakeFiles/pvpre.dir/precon/band.cpp.o
externalpackages/PVODE/libpvpre.a: externalpackages/PVODE/CMakeFiles/pvpre.dir/build.make
externalpackages/PVODE/libpvpre.a: externalpackages/PVODE/CMakeFiles/pvpre.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/ylang/DATA/I-mode/BOUT-dev/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libpvpre.a"
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && $(CMAKE_COMMAND) -P CMakeFiles/pvpre.dir/cmake_clean_target.cmake
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pvpre.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
externalpackages/PVODE/CMakeFiles/pvpre.dir/build: externalpackages/PVODE/libpvpre.a

.PHONY : externalpackages/PVODE/CMakeFiles/pvpre.dir/build

externalpackages/PVODE/CMakeFiles/pvpre.dir/clean:
	cd /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE && $(CMAKE_COMMAND) -P CMakeFiles/pvpre.dir/cmake_clean.cmake
.PHONY : externalpackages/PVODE/CMakeFiles/pvpre.dir/clean

externalpackages/PVODE/CMakeFiles/pvpre.dir/depend:
	cd /media/ylang/DATA/I-mode/BOUT-dev && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/ylang/DATA/I-mode/BOUT-dev /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE /media/ylang/DATA/I-mode/BOUT-dev /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE /media/ylang/DATA/I-mode/BOUT-dev/externalpackages/PVODE/CMakeFiles/pvpre.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : externalpackages/PVODE/CMakeFiles/pvpre.dir/depend

