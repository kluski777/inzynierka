CXX=g++

CXXFLGAS=-std=c++11 -Wall

SOURCES = glowny.cpp
TARGET=result

all: $(TARGET)
	./$(TARGET)	

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)
