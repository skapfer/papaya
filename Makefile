
CXXFLAGS += -Ieinclude -g -ggdb
HEADERS = *.h
SUPPORT = util.o marching.o minkval.o readpgm.o

all: marching_test

clean:
	rm -f marching_test *.o

marching_test: $(SUPPORT)
	$(CXX) -o $@ $(SUPPORT)

