#!/usr/bin/env python3
"""Check that `wild-pbwt -t y` finds the same blocks in a transposed matrix (one line per site) as in the usual one.

  test/check_transposed.py [bin] [trials] [seed]

Random panels with wildcards are written as a matrix and as its transpose (also with CRLF line ends, with empty
lines at the end, and read from the standard input). It also checks that the mistakes a file can have are reported.
The exit code is not 0 if anything differs.
"""
import os, random, shutil, subprocess, sys, tempfile

BIN = sys.argv[1] if len(sys.argv) > 1 else "bin/wild-pbwt"
TRIALS = int(sys.argv[2]) if len(sys.argv) > 2 else 200
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)


def run(args, stdin=None):
    return subprocess.run([BIN] + args, input=stdin, capture_output=True, text=True)


def blocks(args, stdin=None):
    r = run(args, stdin)
    assert r.returncode == 0, (args, r.stderr[-300:])
    out = [l.rsplit(", ", 2) for l in r.stdout.splitlines()]
    return sorted((tuple(sorted(int(x) for x in rows.strip("[],").split(","))), int(i), int(j)) for rows, i, j in out)


tmp = tempfile.mkdtemp()
try:
    mat, tr = os.path.join(tmp, "p.txt"), os.path.join(tmp, "p.tbm")
    for trial in range(TRIALS):
        t, m, n = random.choice([2, 3, 4, 8]), random.randint(2, 9), random.randint(1, 30)
        rate = random.choice([0, 0.05, 0.2, 0.4])
        panel = ["".join("*" if random.random() < rate else str(random.randrange(t)) for _ in range(n)) for _ in range(m)]
        cols = ["".join(panel[r][c] for r in range(m)) for c in range(n)]
        open(mat, "w").write("".join(row + "\n" for row in panel))
        variants = {"lf": "".join(c + "\n" for c in cols), "crlf": "".join(c + "\r\n" for c in cols),
                    "blank lines": "".join(c + "\n" for c in cols) + "\n\n", "no final newline": "\n".join(cols)}
        for flags in ([], ["-w", "y"]):
            want = blocks(["-f", mat, "-a", str(t), "-o", "y"] + flags)
            for kind, text in variants.items():
                open(tr, "w").write(text)
                assert blocks(["-f", tr, "-t", "y", "-a", str(t), "-o", "y"] + flags) == want, (kind, flags, panel)
                assert blocks(["-f", "-", "-t", "y", "-a", str(t), "-o", "y"] + flags, stdin=text) == want, ("stdin " + kind, flags, panel)

    def rejects(text, needle, args=("-a", "2")):
        r = run(["-f", "-", "-t", "y", *args], text)
        assert r.returncode != 0 and needle in r.stderr, (text, needle, r.returncode, r.stderr[-200:])

    rejects("0101\n010\n", "characters, the first one has 4")
    rejects("0101\n01x1\n", "Unexpected character 'x'")
    rejects("0 10\n", "Unexpected character ' '")
    rejects("0120\n", "use -a 3")
    rejects("\n\n", "No columns")
    r = run(["-f", mat + ".vcf", "-t", "y", "-a", "2"])
    assert r.returncode != 0 and "not for a VCF" in r.stderr
    print("ok: %d random panels, transposed as %s, and the mistakes reported" % (TRIALS, ", ".join(list(variants) + ["stdin"])))
finally:
    shutil.rmtree(tmp)
