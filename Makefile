
CXXFLAGS += -Ieinclude -g -ggdb
HEADERS = *.h
SUPPORT = util.o marching.o minkval.o readpgm.o tinyconf.o readpoly.o \
    label.o

all: aspectstudy

clean:
	rm -f marching_test *.o

# trick to recompile everything when the headers change.
ts.headers: $(HEADERS)
	$(MAKE) clean
	touch $@

marching_test: ts.headers $(SUPPORT)
	$(CXX) -o $@ $(SUPPORT)

aspectstudy: ts.headers $(SUPPORT) aspectstudy.o
	$(CXX) -o $@ $(SUPPORT) aspectstudy.o

.PHONY: all clean
