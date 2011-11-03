OBJS=$(shell echo *.cpp | sed s/\\.cpp/.o/g)
CXX=g++
CXXFLAGS=-ggdb -ffast-math -O3
LD=g++
LDFLAGS=-ggdb
BIN=Sparse

default:
	make --no-print-directory -j12 $(BIN)
clean:
	rm -f $(OBJS) $(BIN)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(VERSIONFLAGS)

Sparse: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)