# create a few PGMs to test the pgm reader (thanks to Matthias)
run:	clean createTestPGM
	./createTestPGM
	convert test-16b-ascii.pgm -depth 16 test-16b-bin.pgm
	convert test-16b-ascii.pgm -depth 8 test-8b-bin.pgm
	convert test-16b-ascii.pgm -compress none -depth 8 test-8b-ascii.pgm
	convert test-16b-ascii.pgm -threshold 50%  test-1b-bin.pbm
	convert test-16b-ascii.pgm -compress none -threshold 50%  test-1b-ascii.pbm
createTestPGM:	createTestPGM.cc
	g++ createTestPGM.cc -o createTestPGM -O2 -Wall -pedantic -ansi

clean:
	rm -f createTestPGM *.pgm *.pbm *~
