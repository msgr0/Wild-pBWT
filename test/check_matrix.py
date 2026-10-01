#!/usr/bin/env python3
"""Check that `wild-pbwt` reads a matrix (one haplotype per line) only if it is one, in both ways of reading it
(memory mapped, and with -g one buffered stream per haplotype).

  test/check_matrix.py [bin]

A valid file is read the same with or without the last newline, and every mistake a file can have is reported
with where it is: a character that is not an allele or *, an allele above -a, lines of different lengths, empty
lines at the end, Windows line ends, an empty file. The exit code is not 0 if something is read that should not be.
"""
import os, shutil, subprocess, sys, tempfile

BIN = sys.argv[1] if len(sys.argv) > 1 else "bin/wild-pbwt"
tmp = tempfile.mkdtemp()
try:
    path = os.path.join(tmp, "m.txt")

    def run(text, *args, stream):
        open(path, "w", newline="").write(text)
        return subprocess.run([BIN, "-f", path, "-a", "2", "-o", "y"] + (["-g", "64"] if stream else []) + list(args), capture_output=True, text=True)

    for stream in (False, True):
        way = "-g" if stream else "mmap"
        good = run("0101\n0*10\n1101\n", stream=stream)
        assert good.returncode == 0 and good.stdout, (way, good.stderr)
        assert run("0101\n0*10\n1101", stream=stream).stdout == good.stdout, (way, "the last newline matters")

        for text, needle in [("0101\n01x1\n0000\n", "Unexpected character 'x' in haplotype 1, column 2"),
                             ("0101\n0 10\n", "Unexpected character ' ' in haplotype 1, column 1"),
                             ("0101\n0210\n", "Allele 2 with alphabet size 2, use -a 3 in haplotype 1, column 1"),
                             ("0101\r\n0110\r\n", "end with \\r"),
                             ("0101\n0110\r\n", "Line 2 does not have 4 characters"),
                             ("0101\n010\n0101\n", "Line 2 does not have 4 characters"),
                             ("0101\n01100\n0101\n", "Line 2 does not have 4 characters"),
                             ("0101\n0110\n\n", "Line 3 does not have 4 characters"),
                             ("0101\n0110\n0101\n\n\n", "Line 4 does not have 4 characters"),
                             ("0101\n0110\n01", "Line 3 does not have 4 characters"),
                             ("", "empty")]:
            r = run(text, stream=stream)
            assert r.returncode != 0 and needle in r.stderr, (way, text, needle, r.returncode, r.stderr[-200:])
            if "haplotype" not in needle:  # a wrong character is found while reading, the rest when the file is opened
                assert r.stdout == "", (way, text, "blocks were printed before the mistake")
    print("ok: a valid matrix is read the same in both ways, and 11 kinds of mistake are reported in each")
finally:
    shutil.rmtree(tmp)
