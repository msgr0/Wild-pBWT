.PHONY: all wild-pbwt gen err clean

all: wild-pbwt gen err
wild-pbwt: bin/wild-pbwt
gen: bin/gen
err: bin/err

CXXFLAGS ?= -O2 -march=native
#uncomment next line if SDSL is intalled under user's home directory
#CXXFLAGS ?= -O2 -march=native -I ~/include -L ~/lib

CPPFLAGS += -MMD
LDLIBS += -lsdsl

bin/gen: src/hap_gen.o
bin/err: src/hap_wild.o
bin/wild-pbwt: src/pbwt.o
bin/gen bin/err bin/wild-pbwt:
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(RM) src/*.o src/*.d bin/*

-include src/*.d
