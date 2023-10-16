CXX=g++
CXXFLGAS=-std=c++11
SOURCES = main.cpp
TARGET=result

all: $(TARGET)
	./$(TARGET)	

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)
