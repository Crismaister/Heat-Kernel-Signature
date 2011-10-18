OBJS=$(shell echo *.cpp | sed s/\\.cpp/.o/g)
CXX=g++
CXXFLAGS=-ggdb
LD=g++
LDFLAGS=-ggdb
BIN=Sparse

default:
	make --no-print-directory -j12 $(BIN)
clean:
	rm $(OBJS)  $(BIN)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

Sparse: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)