#!/usr/bin/env python3
"""Check `wild-pbwt -o y` against the definition of maximal perfect haplotype blocks with wildcards.

  test/check_blocks.py [bin] [trials] [seed] [wm|wild]       random small panels, against a brute force
  test/check_blocks.py [bin] --panel FILE -a T [wm|wild]     one panel of any size, against a DFS over consensus patterns

A block is (K, i, j) with |K| >= 2 such that, on every column i..j, the rows of K do not disagree,
they cannot be extended left or right (columns i-1 and j+1 stop them) and no other row fits them.
  wm   (default) Williams & Mumey: rows of K hold exactly ONE allele besides wildcards on every column
       of the block; a column that stops K holds two alleles or none (only wildcards).
  wild (runs the tool with -w y): rows of K hold AT MOST one allele; only two alleles stop K.
The exit code is not 0 if any block is missing, reported twice, or not a block.
"""
import itertools, random, subprocess, sys, tempfile
from collections import Counter

sys.setrecursionlimit(100000)


def masks(panel, t):
    m, n = len(panel), len(panel[0])
    nw = [[sum(1 << r for r in range(m) if panel[r][c] == str(b)) for b in range(t)] for c in range(n)]
    wl = [sum(1 << r for r in range(m) if panel[r][c] == "*") for c in range(n)]
    return nw, wl


def brute(panel, t, strict):
    """all blocks, by enumerating every set of rows (small panels only)"""
    m, n = len(panel), len(panel[0])
    nw, wl = masks(panel, t)
    sets = [K for K in range(1 << m) if bin(K).count("1") >= 2]
    fits = lambda K, c: sum(1 for b in range(t) if K & nw[c][b]) == 1 if strict else sum(1 for b in range(t) if K & nw[c][b]) <= 1
    ok = [{K for K in sets if fits(K, c)} for c in range(n)]
    out = set()
    for i in range(n):
        cur = set(sets)
        for j in range(i, n):
            cur &= ok[j]
            for K in cur:
                if (i == 0 or K not in ok[i - 1]) and (j == n - 1 or K not in ok[j + 1]) \
                   and not any(not K >> r & 1 and (K | 1 << r) in cur for r in range(m)):
                    out.add((K, i, j))
    return out


def dfs(panel, t, strict):
    """all blocks, by extending a consensus pattern to the left from every end column"""
    m, n = len(panel), len(panel[0])
    nw, wl = masks(panel, t)
    alleles = lambda R, c: sum(1 for b in range(t) if R & nw[c][b])
    stops = (lambda R, c: alleles(R, c) != 1) if strict else (lambda R, c: alleles(R, c) >= 2)
    out = set()
    # without -w the set of rows compatible with the consensus is K itself; with -w, K may have lost
    # the rows that held an allele on some column, and then other rows can still join it
    full = lambda K, i, j: not any(all(alleles(K | 1 << r, c) <= 1 for c in range(i, j + 1)) for r in range(m) if not K >> r & 1)

    def go(i, j, R):
        # only alleles some row of R holds (another one would leave rows out of K that fit it);
        # if there are none, K = R goes on through an all-wildcard column, but only without -w
        rows = {R & (nw[i][b] | wl[i]) for b in range(t) if R & nw[i][b]} or (set() if strict else {R})
        for K in rows:
            if bin(K).count("1") < 2 or (strict and any(not K & ~wl[c] for c in range(i, j + 1))):
                continue  # too small, or some column holds only wildcards (and so would any smaller set)
            if (i == 0 or stops(K, i - 1)) and (j == n - 1 or stops(K, j + 1)) and (strict or full(K, i, j)):
                out.add((K, i, j))
            if i > 0:
                go(i - 1, j, K)

    for j in range(n):
        go(j, j, (1 << m) - 1)
    return out


def run(binary, path, t, strict, *extra):
    flags = [] if strict else ["-w", "y"]
    r = subprocess.run([binary, "-f", path, "-a", str(t), "-o", "y", *flags, *extra], capture_output=True, text=True)
    c = subprocess.run([binary, "-f", path, "-a", str(t), "-c", "y", *flags], capture_output=True, text=True)
    g = subprocess.run([binary, "-f", path, "-a", str(t), "-r", "y", *flags, *extra], capture_output=True, text=True)
    assert r.returncode == 0 and c.returncode == 0 and g.returncode == 0, r.stderr[-300:]
    blocks = []
    for line in r.stdout.splitlines():
        rows, i, j = line.rsplit(", ", 2)
        blocks.append((sum(1 << int(x) for x in rows.strip("[],").split(",")), int(i), int(j)))
    for out in (r, c, g):
        assert f"total_blocks_found\t{len(blocks)}\n" in out.stderr, "the count does not match the blocks output"
    # -r y: the same blocks in the same order, the rows sorted, written as maximal ranges
    ranged = g.stdout.splitlines()
    assert len(ranged) == len(blocks), "-r y writes other blocks"
    for line, (mask, i, j) in zip(ranged, blocks):
        rows, ri, rj = line.rsplit(", ", 2)
        runs = [tuple(map(int, item.partition("-")[::2])) if "-" in item else (int(item), int(item)) for item in rows.strip("[]").split(",")]
        assert all(a <= b for a, b in runs) and all(runs[x][1] + 1 < runs[x + 1][0] for x in range(len(runs) - 1)), ("ranges not sorted or not maximal", line)
        assert (sum(1 << x for a, b in runs for x in range(a, b + 1)), int(ri), int(rj)) == (mask, i, j), ("-r y differs from -o y", line)
    return blocks


def compare(truth, got):
    return len(truth - set(got)), len(set(got) - truth), len(got) - len(set(got))


def main():
    args = sys.argv[1:]
    mode = args.pop() if args and args[-1] in ("wm", "wild") else "wm"
    binary = args.pop(0) if args and not args[0].lstrip("-").isdigit() and args[0] != "--panel" else "bin/wild-pbwt"
    strict = mode == "wm"
    if args and args[0] == "--panel":
        path, t = args[1], int(args[args.index("-a") + 1])
        panel = open(path).read().split()
        truth = dfs(panel, t, strict)
        missing, extra, dups = compare(truth, run(binary, path, t, strict))
        print(f"{path} ({mode}): {len(truth)} blocks, missing {missing}, not a block {extra}, duplicated {dups}")
        sys.exit(bool(missing or extra or dups))

    trials = int(args[0]) if args else 2000
    random.seed(int(args[1]) if len(args) > 1 else 1)
    stats, example = Counter(), None
    for _ in range(trials):
        t, m, n = random.choice([2, 3]), random.randint(2, 8), random.randint(1, 8)
        rate = random.choice([0, 0.1, 0.25, 0.4, 0.6])
        panel = ["".join("*" if random.random() < rate else str(random.randrange(t)) for _ in range(n)) for _ in range(m)]
        truth = brute(panel, t, strict)
        assert truth == dfs(panel, t, strict), f"brute force and DFS oracle disagree on {panel}"
        with tempfile.NamedTemporaryFile("w", suffix=".hap") as f:
            f.write("".join(row + "\n" for row in panel))
            f.flush()
            got = run(binary, f.name, t, strict, *(["-g", "64"] if _ % 4 == 0 else []))
        bad = compare(truth, got)
        key = "with wildcards" if any("*" in row for row in panel) else "no wildcards"
        stats[key, "panels"] += 1
        stats[key, "blocks"] += len(truth)
        stats[key, "panels with exact output"] += not any(bad)
        for name, v in zip(("missing", "not a block", "duplicated"), bad):
            stats[key, name] += v
        if any(bad) and (example is None or m * n < len(example[0]) * len(example[0][0])):
            example = (panel, t, bad)
    print(f"mode {mode}: {'Williams & Mumey blocks' if strict else 'blocks may span all-wildcard columns (-w y)'}")
    for (key, what), v in sorted(stats.items()):
        print(f"{key:15} {what:26} {v}")
    if example:
        print("\nsmallest failing panel (alphabet %d), missing/not a block/duplicated = %s:\n  %s" % (example[1], example[2], "\n  ".join(example[0])))
    sys.exit(example is not None)


main()
