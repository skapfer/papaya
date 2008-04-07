
CXXFLAGS += -Ieinclude -g -ggdb
SUPPORT = util.o marching.o minkval.o

all: marching_test

clean:
	rm -f marching_test *.o

marching_test: $(SUPPORT)
	$(CXX) -o $@ $(SUPPORT)

