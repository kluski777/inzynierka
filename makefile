CXX=g++
CXXFLGAS=-std=c++11
SOURCES = main.cpp
TARGET=compiled

all: $(TARGET)
	time ./$(TARGET)	

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)
