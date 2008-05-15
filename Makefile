
CXXFLAGS += -Ieinclude -g -ggdb -Wall
HEADERS = *.h
SUPPORT = util.o marching.o minkval.o readpgm.o tinyconf.o readpoly.o \
    label.o \
    intersect.o \

BINARIES = aspectstudy marching_test

all: aspectstudy testdriver

clean:
	rm -f $(BINARIES) *.o

# trick to recompile everything when the headers change.
ts.headers: $(HEADERS)
	$(MAKE) clean
	touch $@

testdriver: ts.headers $(SUPPORT) driver.o
	$(CXX) -o $@ $(SUPPORT) driver.o

aspectstudy: ts.headers $(SUPPORT) aspectstudy.o
	$(CXX) -o $@ $(SUPPORT) aspectstudy.o

.PHONY: all clean
