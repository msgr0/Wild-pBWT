.PHONY: all wild-pbwt gen err clean

all: wild-pbwt gen err
wild-pbwt: bin/wild-pbwt
gen: bin/gen
err: bin/err

CXXFLAGS ?= -O2 -march=native

CPPFLAGS += -MMD

# VCF/BCF input needs htslib: found with pkg-config, or make HTSLIB=1 with the flags in CPPFLAGS and LDLIBS
# (-I/-L and -lhts), or make HTSLIB=0 to leave it out
HTSLIB ?= $(shell pkg-config --exists htslib 2>/dev/null && echo 1)
ifeq ($(HTSLIB),1)
src/pbwt.o bin/wild-pbwt: CPPFLAGS += -DWITH_HTSLIB $(shell pkg-config --cflags htslib 2>/dev/null)
bin/wild-pbwt: LDLIBS += $(shell pkg-config --libs htslib 2>/dev/null || echo -lhts)
endif

bin/gen: src/hap_gen.o
bin/err: src/hap_wild.o
bin/wild-pbwt: src/pbwt.o
bin/gen bin/err bin/wild-pbwt:
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(RM) src/*.o src/*.d bin/*

-include src/*.d
