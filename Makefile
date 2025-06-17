# Makefile for C++ project with g++ and Intel oneAPI support

# Compiler options
CXX_GCC = g++
CXX_INTEL = icpx

# Default compiler (g++)
CXX = $(CXX_GCC)

# Directories
SRC_DIR = src
INCLUDE_DIR = include
EXAMPLE_DIR = example
BUILD_DIR = build

# Flags
CXXFLAGS = -Wall -O2 -I$(INCLUDE_DIR)
LDFLAGS =

# Source files and target
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
TARGET = main

# Default target
all: directories $(TARGET)

# Create build directory
directories:
	mkdir -p $(BUILD_DIR)

# Linking
$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

# Compiling
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Intel oneAPI target
intel: CXX = $(CXX_INTEL)
intel: clean all

# Example building (optional)
examples: all
	$(CXX) $(CXXFLAGS) $(EXAMPLE_DIR)/*.cpp -o example_app $(LDFLAGS) -L. -l:$(TARGET)

.PHONY: all clean intel directories examples
