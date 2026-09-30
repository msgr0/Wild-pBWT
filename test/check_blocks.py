#!/usr/bin/env python3
"""Compare `wild-pbwt -o y` with a brute-force MPHBw enumeration on random small panels.

usage: test/check_blocks.py [path/to/wild-pbwt] [trials] [seed]

MPHBw (Williams & Mumey): (K, i, j), |K| >= 2, such that on columns i..j the rows of K are
pairwise compatible (equal or '*'), and the block cannot be extended left, right, or by a row.
"""
import itertools, random, subprocess, sys, tempfile
from collections import Counter

BIN = sys.argv[1] if len(sys.argv) > 1 else "bin/wild-pbwt"
TRIALS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)


def compatible(panel, K, c):
    return len({panel[r][c] for r in K} - {"*"}) <= 1


def classify(panel, K, i, j):
    m, n = len(panel), len(panel[0])
    if not all(compatible(panel, K, c) for c in range(i, j + 1)):
        return "invalid"
    why = []
    if i > 0 and compatible(panel, K, i - 1):
        why.append("not left-maximal")
    if j < n - 1 and compatible(panel, K, j + 1):
        why.append("not right-maximal")
    if any(all(compatible(panel, K | {r}, c) for c in range(i, j + 1)) for r in set(range(m)) - K):
        why.append("not row-maximal")
    return " + ".join(why) or "ok"


def brute(panel):
    m, n = len(panel), len(panel[0])
    return {(frozenset(K), i, j)
            for size in range(2, m + 1) for K in itertools.combinations(range(m), size)
            for i in range(n) for j in range(i, n)
            if classify(panel, frozenset(K), i, j) == "ok"}


def tool(panel, t):
    with tempfile.NamedTemporaryFile("w", suffix=".hap") as f:
        f.write("".join(row + "\n" for row in panel))
        f.flush()
        r = subprocess.run([BIN, "-f", f.name, "-a", str(t), "-o", "y"], capture_output=True, text=True)
    assert r.returncode == 0, (panel, r.stderr)
    blocks = []
    for line in r.stdout.splitlines():
        rows, i, j = line.rsplit(", ", 2)
        blocks.append((frozenset(int(x) for x in rows.strip("[],").split(",")), int(i), int(j)))
    return blocks


stats, examples = Counter(), {}
for _ in range(TRIALS):
    t, m, n = random.choice([2, 3]), random.randint(2, 6), random.randint(1, 7)
    rate = random.choice([0, 0.1, 0.25, 0.4])
    panel = ["".join("*" if random.random() < rate else str(random.randrange(t)) for _ in range(n)) for _ in range(m)]
    truth, got = brute(panel), tool(panel, t)
    key = "no wildcards" if not any("*" in row for row in panel) else "with wildcards"
    issues = ["missing"] * len(truth - set(got)) + ["duplicate"] * (len(got) - len(set(got)))
    issues += [why for why in (classify(panel, *b) for b in set(got)) if why != "ok"]
    stats[key, "panels"] += 1
    stats[key, "true blocks"] += len(truth)
    stats[key, "reported lines"] += len(got)
    stats[key, "panels with exact output"] += not issues
    for why in issues:
        stats[key, why] += 1
        if why not in examples or m * n < len(examples[why][0]) * len(examples[why][0][0]):
            examples[why] = (panel, t, truth, got)

for (key, what), v in sorted(stats.items()):
    print(f"{key:15} {what:45} {v}")
for why, (panel, t, truth, got) in examples.items():
    fmt = lambda bs: sorted((sorted(K), i, j) for K, i, j in bs)
    print(f"\nsmallest example of '{why}' (alphabet {t}):\n  " + "\n  ".join(panel))
    print("  true:    ", fmt(truth), "\n  reported:", fmt(got))
sys.exit(any(what in examples for _, what in stats))
