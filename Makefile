all: pbwt
debug: debug-pbwt
#-DNDEBUG
SDSL:= -I ~/include -L ~/lib
SDSLPOST:= -lsdsl
FLAGS:= -Ofast $(SDSL)
DEFLAGS:= -DDEBUG -Ofast -Wall
FLAGSHAPLO:= -Ofast -std=c++11

old: src/pbwt-old.cpp src/MatrixReader.hpp
	c++ $(FLAGS) src/pbwt-old.cpp -o bin/old $(SDSLPOST)
new: src/pbwt-new.cpp src/MatrixReader.hpp
	c++ $(FLAGS) src/pbwt-new.cpp -o bin/new $(SDSLPOST)
gen: src/hap_gen.cpp
	c++ $(FLAGS) src/hap_gen.cpp -o bin/gen 
err: src/hap_wild.cpp
	c++ $(FLAGS) src/hap_wild.cpp -o bin/err
pbwt: src/pbwt.cpp src/FileReader.hpp
	c++ $(FLAGS) src/pbwt.cpp -o bin/pbwt $(SDSLPOST)
debug-pbwt: src/pbwt.cpp src/MatrixReader.hpp
	c++ $(DEFLAGS) src/pbwt.cpp -o bin/pbwt
clean:
	rm bin/*
hap: src/haplotype-pbwt-lite.cpp
	g++ $(FLAGSHAPLO) src/haplotype-pbwt-lite.cpp -o bin/haplotype-pbwt-lite



