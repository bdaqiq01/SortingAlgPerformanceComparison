# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 $(shell python3-config --cflags)
LDFLAGS = $(shell python3-config --ldflags) -lpython3.10

# Output executable name
TARGET = main

# Source files
SRCS = main.cpp

# Header files
HEADERS = plotSort.h matplotlibcpp.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule to build the executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) -lm

# Rule for object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJS) $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)
