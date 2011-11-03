OBJS=$(shell echo *.cpp | sed s/\\.cpp/.o/g)
CXX=g++
NUTRLAN_DIR=nutrlan-0.2
CXXFLAGS=-ggdb -ffast-math -O3 -DBLAS
LD=g++
LDFLAGS=-ggdb -L$(NUTRLAN_DIR) -lnutrlan -llapack -lblas
BIN=Sparse

INCFLAGS=-I$(NUTRLAN_DIR)/INC

default:
	make --no-print-directory -j12 $(BIN)
clean:
	rm -f $(OBJS) $(BIN)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@ $(VERSIONFLAGS)

Sparse: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)