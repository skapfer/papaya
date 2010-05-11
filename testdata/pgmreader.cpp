// vim: et:sw=4:ts=4

#include <iostream>
#include "../util.h"

int main (int argc, char **argv) {

    bool failed = false;
    Pixmap ascii, binary;

    std::cerr << "Testing PGM reader...\n";

    load_pgm (&ascii, "pgmreader_test_in/odd/test-16b-ascii.pgm");
    load_pgm (&binary, "pgmreader_test_in/odd/test-16b-bin.pgm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/odd/test-16b: ascii / binary mismatch\n";
        failed = true;
    }

    load_pgm (&ascii, "pgmreader_test_in/odd/test-8b-ascii.pgm");
    load_pgm (&binary, "pgmreader_test_in/odd/test-8b-bin.pgm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/odd/test-8b: ascii / binary mismatch\n";
        failed = true;
    }

    load_pgm (&ascii, "pgmreader_test_in/odd/test-1b-ascii.pbm");
    load_pgm (&binary, "pgmreader_test_in/odd/test-1b-bin.pbm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/odd/test-1b: ascii / binary mismatch\n";
        failed = true;
    }

    load_pgm (&ascii, "pgmreader_test_in/even/test-16b-ascii.pgm");
    load_pgm (&binary, "pgmreader_test_in/even/test-16b-bin.pgm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/even/test-16b: ascii / binary mismatch\n";
        failed = true;
    }

    load_pgm (&ascii, "pgmreader_test_in/even/test-8b-ascii.pgm");
    load_pgm (&binary, "pgmreader_test_in/even/test-8b-bin.pgm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/even/test-8b: ascii / binary mismatch\n";
        failed = true;
    }

    load_pgm (&ascii, "pgmreader_test_in/even/test-1b-ascii.pbm");
    load_pgm (&binary, "pgmreader_test_in/even/test-1b-bin.pbm");
    if (! (ascii == binary)) {
        std::cerr << "pgmreader_test_in/even/test-1b: ascii / binary mismatch\n";
        failed = true;
    }

    if (!failed) std::cerr << "OK\n";

    return int (failed);
}
