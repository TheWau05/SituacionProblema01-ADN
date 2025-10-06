# Compiler
CXX = g++

# Compiler flags
# -std=c++17: Use the C++17 standard
# -Wall: Enable all warnings
# -Wextra: Enable extra warnings
# -g: Generate debugging information
CXXFLAGS = -std=c++17 -Wall -Wextra -g

# The name of the final executable
TARGET = main_executable

# List of all source (.cpp) files
SOURCES = SituacionP1.cpp ReadFasta.cpp LongestPalindrome.cpp

# Automatically generate object (.o) file names from source files
OBJECTS = $(SOURCES:.cpp=.o)

# Default rule: build the target executable
all: $(TARGET)

# Rule to link the object files into the final executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Rule to compile a .cpp file into a .o file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up generated files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets are not actual files
.PHONY: all clean
