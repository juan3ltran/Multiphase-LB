CXX = g++
CXXFLAGS = -std=c++11 -fopenmp -Wall -I. -O3
LDFLAGS = -fopenmp
SOURCES = main.cpp Latticeboltzmann.cpp animation.cpp vector.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = Multiphase
GIF = phases.gif

all: $(GIF)

$(GIF):$(EXECUTABLE)
	./$< |gnuplot

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ 


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.dat 