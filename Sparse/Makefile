OBJS=$(shell echo *.cpp | sed s/\\.cpp/.o/g)
CXX=g++
NUTRLAN_DIR=nutrlan-0.2
NUTRLAN_TAR=nutrlan-0.2.tar.gz
CXXFLAGS=-ggdb -ffast-math -O3 -DBLAS
LD=g++
LDFLAGS=-ggdb -L$(NUTRLAN_DIR) -L. -lnutrlan -llapack -lblas
BIN=Sparse

INCFLAGS=-I$(NUTRLAN_DIR)/INC

default:
	make --no-print-directory -j12 $(BIN)
clean:
	rm -f $(OBJS) $(BIN) $(NUTRLAN_DIR)/Makefile

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@ $(VERSIONFLAGS)

Sparse: $(OBJS) libnutrlan.a
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

main.o: main.cpp $(NUTRLAN_DIR)/Makefile

$(NUTRLAN_DIR)/Makefile:
	tar --exclude Make.inc -zxf $(NUTRLAN_DIR)/$(NUTRLAN_TAR)

libnutrlan.a: $(NUTRLAN_DIR)/Makefile
	+make -C $(NUTRLAN_DIR)
