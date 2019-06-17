
LIBRARIES =  -lblas -llapack
CXX = g++
CXXFLAGS = -O3


optimalFlow:	main.cpp vector.h foundfeasible.h foundfeasible.o checkedcolumnstree.o
	$(CXX) $(CXXFLAGS) -o $@ main.cpp foundfeasible.o checkedcolumnstree.o $(LIBRARIES)

clean:
	rm -fr *.o optimalFlow
