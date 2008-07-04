
CXXFLAGS += -Ieinclude -g -ggdb -Wall
HEADERS = *.h
SUPPORT = util.o marching.o minkval.o readpgm.o tinyconf.o readpoly.o \
    label.o \
    intersect.o \

BINARIES = papaya testdata/eigensystem testdata/tsvdiff

all: $(BINARIES)

clean:
	rm -f $(BINARIES) *.o

# trick to recompile everything when the headers change.
ts.headers: $(HEADERS)
	$(MAKE) clean
	touch $@

papaya: ts.headers $(SUPPORT) driver.o
	$(CXX) -o $@ $(SUPPORT) driver.o

testdata/tsvdiff: ts.headers util.o tsvdiff.o
	$(CXX) -o $@ util.o tsvdiff.o

testdata/eigensystem: ts.headers $(SUPPORT) testdata/eigensystem.cpp
	$(CXX) $(CXXFLAGS) -o $@ $(SUPPORT) testdata/eigensystem.cpp

.PHONY: all clean
