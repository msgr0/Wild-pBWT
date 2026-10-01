#!/usr/bin/env python3
"""Check that `wild-pbwt` finds the same blocks in a VCF/BCF as in the same panel written as a matrix.

  test/check_vcf.py [bin] [trials] [seed]

Random panels with wildcards are written as a matrix, as a VCF (also read from the standard input) and,
if bgzip and bcftools are installed, as .vcf.gz and .bcf. Haplotype h of sample s is the row ploidy*s+h of the matrix.
It also checks that the mistakes a VCF can have are reported. The exit code is not 0 if anything differs.
"""
import os, random, shutil, subprocess, sys, tempfile

BIN = sys.argv[1] if len(sys.argv) > 1 else "bin/wild-pbwt"
TRIALS = int(sys.argv[2]) if len(sys.argv) > 2 else 60
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
ALTS = ["C", "G", "T", "CC", "GG", "TT", "AC"]


def vcf_text(panel, t, ploidy=2, phased=True, extra=False, format_tags="GT"):
    names = ["s%d" % i for i in range(len(panel) // ploidy)]
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=1>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)]
    for c in range(len(panel[0])):
        cells = [("|" if phased else "/").join("." if panel[s * ploidy + h][c] == "*" else panel[s * ploidy + h][c] for h in range(ploidy))
                 + (":7" if extra else "") for s in range(len(names))]
        if format_tags != "GT":
            cells = ["7"] * len(names)
        fmt = format_tags if format_tags != "GT" else ("GT:DP" if extra else "GT")
        lines.append("\t".join(["1", str(100 + 10 * c), ".", "A", ",".join(ALTS[:t - 1]), ".", ".", ".", fmt] + cells))
    return "\n".join(lines) + "\n"


def run(args, stdin=None):
    return subprocess.run([BIN] + args, input=stdin, capture_output=True, text=True)


def blocks(args, stdin=None):
    r = run(args, stdin)
    assert r.returncode == 0, (args, r.stderr[-300:])
    out = [l.rsplit(", ", 2) for l in r.stdout.splitlines()]
    return sorted((tuple(sorted(int(x) for x in rows.strip("[],").split(","))), int(i), int(j)) for rows, i, j in out), r.stderr


if "cannot read VCF" in run(["-f", "x.vcf", "-a", "2"]).stderr:
    print("this build has no htslib: nothing to check")
    sys.exit(0)
tmp = tempfile.mkdtemp()
try:
    mat, vcf = os.path.join(tmp, "p.txt"), os.path.join(tmp, "p.vcf")
    for trial in range(TRIALS):
        t, samples, n = random.choice([2, 3, 4]), random.randint(1, 6), random.randint(1, 30)
        rate = random.choice([0, 0.05, 0.2, 0.4])
        panel = ["".join("*" if random.random() < rate else str(random.randrange(t)) for _ in range(n)) for _ in range(2 * samples)]
        open(mat, "w").write("".join(row + "\n" for row in panel))
        text = vcf_text(panel, t, extra=trial % 2 == 0)
        open(vcf, "w").write(text)
        files = {"vcf": vcf}
        if shutil.which("bgzip"):
            subprocess.run("bgzip -c %s > %s.gz" % (vcf, vcf), shell=True, check=True)
            files["vcf.gz"] = vcf + ".gz"
        if shutil.which("bcftools"):
            subprocess.run("bcftools view -O b -o %s.bcf %s" % (vcf[:-4], vcf), shell=True, check=True)
            files["bcf"] = vcf[:-4] + ".bcf"
        for flags in ([], ["-w", "y"]):
            want, _ = blocks(["-f", mat, "-a", str(t), "-o", "y"] + flags)
            for kind, path in files.items():
                assert blocks(["-f", path, "-a", str(t), "-o", "y"] + flags)[0] == want, (kind, flags, panel)
            assert blocks(["-f", "-", "-a", str(t), "-o", "y"] + flags, stdin=text)[0] == want, ("stdin", flags, panel)
        if t == 2:  # without -a a VCF has two alleles
            assert blocks(["-f", vcf, "-o", "y"])[0] == blocks(["-f", mat, "-a", "2", "-o", "y"])[0], panel

    hap = ["0101*1", "01*011", "1101*0"]  # haploid samples: one row each
    open(mat, "w").write("\n".join(hap) + "\n")
    assert blocks(["-f", "-", "-a", "2", "-o", "y"], stdin=vcf_text(hap, 2, ploidy=1))[0] == blocks(["-f", mat, "-a", "2", "-o", "y"])[0]

    panel = ["012", "102", "0*1", "110"]
    r = run(["-f", "-", "-a", "2"], vcf_text(panel, 3))
    assert r.returncode != 0 and "use -a 3" in r.stderr, "an allele above the alphabet size is not reported"
    r = run(["-f", "-", "-a", "2"], vcf_text(panel, 2, format_tags="DP"))
    assert r.returncode != 0 and "No genotypes" in r.stderr, "a record without genotypes is not reported"
    r = run(["-f", "-", "-a", "2", "-c", "y"], vcf_text(["01", "10", "01", "01"], 2, phased=False))
    assert r.returncode == 0 and "unphased" in r.stderr, "unphased genotypes are not reported"
    print("ok: %d random panels as %s, stdin, haploid, and the mistakes reported" % (TRIALS, ", ".join(files)))
finally:
    shutil.rmtree(tmp)
