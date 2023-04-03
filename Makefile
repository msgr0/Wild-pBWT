.PHONY: all wild-pbwt gen err clean depend

all: wild-pbwt gen err
wild-pbwt: bin/wild-pbwt
gen: bin/gen
err: bin/err

CXXFLAGS ?= -O2 -march=native 
#uncomment next line if SDSL is intalled under user's home directory
#CXXFLAGS ?= -O2 -march=native -I ~/include -L ~/lib


LDLIBS += -lsdsl


SRCS=src/hap_gen.cpp src/hap_wild.cpp src/pbwt.cpp

bin/gen: src/hap_gen.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

bin/err: src/hap_wild.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

bin/wild-pbwt: src/pbwt.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

bin/debug-wild-pbwt: src/pbwt.cpp
	$(CXX) $(DEBUG_CXXFLAGS) $(DEBUG_CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(RM) src/*.o bin/* .depend

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^ >>./.depend;

include .depend
