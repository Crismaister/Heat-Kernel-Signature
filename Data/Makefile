
FILES=$(shell echo *.[oO][fF][fF])
SFILES=$(shell echo *.[oO][fF][fF] | sed 's/.[oO][fF][fF]/.off.signature.$(VERSION)/g')

default: $(SFILES)

%.off.signature.$(VERSION) : %.off
	../Sparse/Sparse $< $(VERSION)
#	mv $<.signature $<.signature.$(VERSION)
