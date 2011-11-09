CXX = g++
CXXFLAGS := -c
INCLUDES := -I.	-I/home/lfdguest000/oreccb/trunk

LDFLAGS  :=
OBJECTS  := benchmark.o MyPriorityQueue.o presentation1.o presentation2.o runSim.o sortfile.o main.o
EXECUTABLE :=tbiexe


$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@

clean:
	rm *.o tbiexe 
