#include "MatrixReader.hpp"
#include "TransposedReader.hpp"
#ifdef WITH_HTSLIB
#include "VcfReader.hpp"
#endif
#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

typedef long long int int_t;

// globals
bool verbose = false;
bool count_blocks = false;
bool output_blocks = false;
bool wild_columns = false; // -w y: blocks may span (and stop at) columns where all their rows have a wildcard
int_t minimal_block_size = 2;

class PbwtOrder
{
private:
    std::vector<int_t> ak;
    std::vector<int_t> dk;
    struct OpenNode
    {
        int_t d, m;
        size_t ties; // where its ties start in the stack ties
    };
    std::vector<OpenNode> open_nodes; // see gap_blocks
    std::vector<int_t> ties;
    std::vector<std::vector<int_t>> count; // count[l][x]: entries ak[0..x-1] with allele l on the next column
    std::vector<std::vector<uint64_t>> wild_row; // per row, one bit per column so far: is it a wildcard
    std::vector<int_t> block_rows;
    std::vector<int_t> in_block; // per row, last block/group it was seen in (see maximal_rows)
    std::vector<int_t> in_group;
    int_t stamp = 0;
    std::vector<std::vector<int_t>> where; // per row, the positions of its copies in ak

    int_t k = 0;
    bool more = false; // there is a column after the current one
    int_t M;
    const allele_t *column_pointer;

    allele_t alphabet_size;

public:
    int_t total_blocks = 0;
    int_t expansion_count = 0;
    int_t collapsed_rows_count = 0;
    int_t collapse_count = 0;

    PbwtOrder(const allele_t *r, int_t lines, allele_t alphabet_size)
        : ak(lines), dk(lines, 0), wild_row(lines), in_block(lines, 0), in_group(lines, 0), where(lines), M(lines), column_pointer(r), alphabet_size(alphabet_size)
    {
        std::iota(ak.begin(), ak.end(), 0);
    }
    void see()
    {
        auto cell = [&](int_t i)
        { return column_pointer[i] == WILD ? "*" : std::to_string(column_pointer[i]); };

        std::cerr << "k: " << k << "\ncolumn: \t";
        for (int_t i = 0; i < M; i++)
            std::cerr << cell(i) << " ";
        std::cerr << "\nordered: \t";
        for (auto &i : ak)
            std::cerr << cell(i) << " ";
        std::cerr << "\nak: \t\t";
        for (auto &i : ak)
            std::cerr << i << " ";
        std::cerr << "\ndk: \t\t";
        for (auto &i : dk)
            std::cerr << i << " ";
        std::cerr << "\n------------------------\n";
    }

    int_t get_current_size() { return ak.size(); }

    bool wild_at(int_t row, int_t column) { return wild_row[row][column >> 6] >> (column & 63) & 1; }

    // How many different alleles (wildcards excluded) the entries ak[m..n] have on column k + 1:
    // 0, 1, or 2 for two or more. Like in the original haploblocks code the answer comes from prefix
    // sums over the current order (one per allele), so it costs O(alleles) for any interval.
    int_t next_alleles(int_t m, int_t n)
    {
        int_t seen = 0;
        for (int_t l = 0; l < alphabet_size && seen < 2; l++)
        {
            seen += count[l][n + 1] > count[l][m];
        }
        return seen;
    }

    // How many different alleles (wildcards excluded), 2 for two or more, the rows of the entries
    // ak[m..n] have on column d - 1. The entries split where dk == d (the positions in ties) into
    // groups that agree on that column: a row with an allele there is in one group only, a row with
    // a wildcard is in all of them. So it is enough to find in each group a row that is not a
    // wildcard, and that is nearly always its first entry.
    int_t previous_alleles(int_t m, int_t n, int_t d, const int_t *ties, int_t ties_count)
    {
        int_t alleles = 0;
        for (int_t g = 0; g <= ties_count && alleles < 2; g++)
        {
            int_t from = g == 0 ? m : ties[g - 1];
            int_t to = g == ties_count ? n : ties[g] - 1;
            for (int_t a = from; a <= to; a++)
            {
                if (!wild_at(ak[a], d - 1))
                {
                    alleles++;
                    break;
                }
            }
        }
        return alleles;
    }

    // Williams & Mumey blocks (K, d, k): on every column d..k the rows of K hold exactly one
    // allele besides wildcards (so no column with only wildcards inside the block), and at
    // columns d - 1 and k + 1 they do not (two alleles, or only wildcards, both stop the block).
    // The rows of K are the distinct rows of an interval ak[m..n] of entries matching from d, and
    // that interval holds every row compatible with the block, so K is row-maximal by construction.
    // The right end has been checked by the caller.
    bool williams_mumey_block(int_t m, int_t n, int_t d, const int_t *ties, int_t ties_count)
    {
        if (d > 0 && previous_alleles(m, n, d, ties, ties_count) == 1)
            return false;
        for (int_t c = d; c <= k; c++) // every column needs a row with an allele
        {
            int_t a = m;
            while (a <= n && wild_at(ak[a], c))
                a++;
            if (a > n)
                return false;
        }
        return true;
    }

    // The distinct rows of ak[m..n], in block_rows.
    int_t collect_rows(int_t m, int_t n)
    {
        int_t block = ++stamp;
        block_rows.clear();
        for (int_t a = m; a <= n; a++)
        {
            if (in_block[ak[a]] != block)
            {
                in_block[ak[a]] = block;
                block_rows.push_back(ak[a]);
            }
        }
        return block_rows.size();
    }

    // With -w y, where a column with only wildcards does not stop a block: the number of distinct
    // rows of ak[m..n] (left in block_rows), or 0 if the block (rows, d, k) must not be reported
    // because it is not row-maximal or it is reported from elsewhere. The left end has been checked.
    //
    // When all the rows have a wildcard on the same column, every other filling of it gives another
    // run of entries matching from d, with the same rows plus the ones that really have that
    // allele. Any such run holds a copy of each row of the block, so the copies of one row find
    // them all: a run with more rows means the block is not row-maximal, a run with the same rows
    // is the same block, reported from its first run.
    int_t maximal_rows(int_t m, int_t n, int_t d)
    {
        int_t range = collect_rows(m, n), block = stamp;
        int_t fewest = ak[m]; // row of the block with the fewest copies
        for (int_t r : block_rows)
        {
            if (where[r].size() < where[fewest].size())
                fewest = r;
        }

        int_t hi = 0;
        for (int_t p : where[fewest])
        {
            if ((p >= m && p <= n) || p < hi)
                continue;
            int_t lo = p;
            while (lo > 0 && dk[lo] <= d)
                lo--;
            for (hi = p + 1; hi < get_current_size() && dk[hi] <= d; hi++)
                ;
            int_t run = ++stamp, rows = 0, shared = 0;
            for (int_t a = lo; a < hi; a++)
            {
                if (in_group[ak[a]] != run)
                {
                    in_group[ak[a]] = run;
                    rows++;
                    shared += in_block[ak[a]] == block;
                }
            }
            if (shared == range && (rows > range || lo < m))
                return 0;
        }
        return range;
    }

    // A block found on the entries ak[m..n], which match from column d to column k, and where
    // the ties are the positions inside (m, n] with dk == d.
    void report_block(int_t m, int_t n, int_t d, const int_t *ties, int_t ties_count)
    {
        if (n - m == 1 && ak[m] == ak[n]) // one row (a row is not twice in a row in ak more than twice)
            return;
        if (more) // right: one allele on column k + 1 means the block goes on, and so do
        {              // no alleles unless -w y (two alleles is what stops it there)
            int_t alleles = next_alleles(m, n);
            if (wild_columns ? alleles < 2 : alleles == 1)
                return;
        }
        int_t range = 0;
        if (wild_columns)
        {
            if (d > 0 && previous_alleles(m, n, d, ties, ties_count) < 2)
                return;
            range = maximal_rows(m, n, d);
            if (range == 0)
                return;
        }
        else
        {
            if (!williams_mumey_block(m, n, d, ties, ties_count))
                return;
            if (!count_blocks)
                range = collect_rows(m, n);
        }
        if (count_blocks)
        {
            total_blocks += 1;
        }
        else if ((k - d + 1) * range >= minimal_block_size)
        {
            total_blocks++;
            if (output_blocks)
            {
                std::cout << "[";
                for (int_t r : block_rows)
                {
                    std::cout << r << ",";
                }
                std::cout << "], " << d << ", " << k << "\n";
            }
        }
    }

    // The blocks ending at column k start where dk says so: each is an interval ak[m..n] whose
    // entries match from column d, with d the largest dk inside (m, n] and dk[m], dk[n + 1] larger.
    // Like update_interval in haploblocks, one pass with a stack finds them all: the stack holds
    // the intervals still open, with decreasing d, each growing while dk stays at most its d. The
    // positions where dk equals d are kept too (in ties, also a stack): they split the interval
    // in the groups that differ on column d - 1.
    // Reports the blocks that end at the current column. If more_columns, the column after it is
    // in column_pointer, what tells if a block can be extended to the right.
    void gap_blocks(bool more_columns)
    {
        more = more_columns;
        int_t current_size = get_current_size();
        if (more)
        {
            count.assign(alphabet_size, std::vector<int_t>(current_size + 1, 0));
            for (int_t x = 0; x < current_size; x++)
            {
                for (int_t l = 0; l < alphabet_size; l++)
                {
                    count[l][x + 1] = count[l][x];
                }
                allele_t v = column_pointer[ak[x]];
                if (v != WILD)
                {
                    count[v][x + 1]++;
                }
            }
        }
        if (wild_columns)
        {
            for (auto &w : where)
            {
                w.clear();
            }
            for (int_t x = 0; x < current_size; x++)
            {
                where[ak[x]].push_back(x);
            }
        }
        open_nodes.clear();
        ties.clear();
        for (int_t i = 1; i <= current_size; i++)
        {
            int_t v = i < current_size ? dk[i] : std::numeric_limits<int_t>::max(); // nothing matches past the end
            int_t m = i - 1;
            while (!open_nodes.empty() && open_nodes.back().d < v)
            {
                const OpenNode &o = open_nodes.back();
                if (o.d <= k) // a larger d means the entries differ on column k
                {
                    report_block(o.m, i - 1, o.d, ties.data() + o.ties, ties.size() - o.ties);
                }
                m = o.m;
                ties.resize(o.ties);
                open_nodes.pop_back();
            }
            if (i < current_size)
            {
                if (open_nodes.empty() || open_nodes.back().d != v)
                {
                    open_nodes.push_back({v, m, ties.size()});
                }
                ties.push_back(i);
            }
        }
        this->k++;
    }

    void next()
    {
        int_t size_curr = 0;
        std::vector<std::vector<int_t>> a(alphabet_size);
        std::vector<std::vector<int_t>> d(alphabet_size);
        std::vector<int_t> p(alphabet_size, k + 1);

        for (int_t r = 0; r < M; r++)
        {
            if ((k & 63) == 0)
            {
                wild_row[r].push_back(0);
            }
            if (column_pointer[r] == WILD)
            {
                wild_row[r][k >> 6] |= uint64_t(1) << (k & 63);
            }
        }

        for (int_t l = 0; l < alphabet_size; l++)
        {
            a[l].reserve(get_current_size());
            d[l].reserve(get_current_size());
        }

        auto push = [&](int_t l, int_t i)
        {
            a[l].emplace_back(ak[i]);
            d[l].emplace_back(p[l]);
            p[l] = 0;
            size_curr += 1;
        };

        for (int_t i = 0; i < get_current_size(); i++)
        {
            allele_t allele = column_pointer[ak[i]];

            for (int_t l = 0; l < alphabet_size; l++)
            {
                if (dk[i] > p[l])
                {
                    p[l] = dk[i];
                }
            }
            if (allele != WILD)
            {
                push(allele, i);
            }
            else
            {
                this->expansion_count += 1;
                for (int_t l = 0; l < alphabet_size; l++)
                {
                    push(l, i);
                }
            }
        }

        ak.clear();
        dk.clear();
        ak.reserve(size_curr);
        dk.reserve(size_curr);

        for (int_t i = 0; i < alphabet_size; i++)
        {
            ak.insert(ak.end(), a[i].begin(), a[i].end());
            dk.insert(dk.end(), d[i].begin(), d[i].end());
        }

        std::vector<int_t> ac;
        std::vector<int_t> dc;
        ac.reserve(size_curr);
        dc.reserve(size_curr);

        bool started = false;
        int_t max_d = 0;
        int_t upper_d = 0;
        for (int_t i = 0; i < size_curr - 1; i++)
        {
            if (ak[i] == ak[i + 1])
            {
                if (!started)
                {
                    started = true;
                    upper_d = dk[i];
                    ac.emplace_back(ak[i]);
                    dc.emplace_back(upper_d);
                }
                else
                {
                    this->collapsed_rows_count += 1;
                }
                max_d = std::max(dk[i + 1], max_d);
            }
            else
            {
                if (started)
                {
                    // The run is the copies b..i of one row, with upper_d the divergence of the first
                    // one from the entry before the run, dk[i + 1] the one of the entry after it from the
                    // last copy, and max_d the largest between two copies. Entries of the same row never
                    // get separated later, and divergences only grow, so one entry is enough unless some
                    // d has max_d > d >= both ends: then the entry before and the one after would
                    // wrongly be put in one block through this row. Only then the last copy is kept.
                    if (max_d <= std::max(upper_d, dk[i + 1]))
                    {
                        this->collapsed_rows_count += 1;
                    }
                    else
                    {
                        ac.emplace_back(ak[i]);
                        dc.emplace_back(max_d);
                    }
                    this->collapse_count += 1;
                    max_d = 0;
                    started = false;
                }
                else
                {
                    ac.emplace_back(ak[i]);
                    dc.emplace_back(dk[i]);
                }
            }
        }

        if (!started)
        {
            ac.emplace_back(ak[size_curr - 1]);
            dc.emplace_back(dk[size_curr - 1]);
        }
        else
        {
            this->collapsed_rows_count += 1;
        }

        ak = std::move(ac);
        dk = std::move(dc);

        if (verbose)
        {
            see();
        }
    }
};

void usage()
{
    std::cerr << "Usage: ./bin/wild-pbwt \n";
    std::cerr << "-a <alphabet_size> \n-f <filename> <a matrix, or a .vcf .vcf.gz .bcf file, or - for a VCF/BCF on stdin> \n";
    std::cerr << "-t y <the matrix is transposed: one line per site, one character per haplotype; - is stdin> \n";
    std::cerr << "-c y <count max blocks> \n-o y <out_blocks to std_out> \n";
    std::cerr << "-b <min block size> -v y <verbose> \n";
    std::cerr << "-g <buffer_size> <read with one buffered stream per row instead of mmap> \n";
    std::cerr << "-w y <blocks may span columns where all their rows have a wildcard, each block once>\n";
    std::cerr << "     default: Williams & Mumey blocks, exactly one allele per column of a block \n";
    std::cerr << "Example: ./bin/wild-pbwt -a 3 -f paper_wild -o y > paper_wild_out.txt \n";
    exit(0);
}

// A VCF or BCF, by its name, or the standard input, is read with htslib; anything else is a matrix.
bool is_vcf(const std::string &name)
{
    for (const char *ext : {".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"})
    {
        size_t n = strlen(ext);
        if (name.size() >= n && name.compare(name.size() - n, n, ext) == 0)
            return true;
    }
    return name == "-";
}

int main(int argc, char **argv)
{
    int ch;
    std::string filename;
    allele_t alphabet_size = 0;
    int_t buffer_size = 0;
    bool transposed = false;

    while ((ch = getopt(argc, argv, "hf:a:c:o:v:w:t:g:b:")) != -1)
    {
        switch (ch)
        {
        case 'h':
            usage();
            break;
        case 'f':
            filename = optarg;
            break;
        case 'a':
            alphabet_size = atoi(optarg);
            break;
        case 'c':
            count_blocks = true;
            break;
        case 'o':
            output_blocks = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 'w':
            wild_columns = true;
            break;
        case 't':
            transposed = true;
            break;
        case 'g':
            buffer_size = atoi(optarg);
            if (buffer_size < 2)
            {
                std::cerr << "buffersize is too small\n";
                exit(0);
            }
            break;
        case 'b':
            minimal_block_size = atoi(optarg);
            break;
        case '?':
        default:
            usage();
        }
    }
    if (filename.empty())
        usage();
    bool vcf = is_vcf(filename);
    if (transposed && vcf && filename != "-")
    {
        std::cerr << "-t is for a text matrix, not for a VCF/BCF\n";
        exit(1);
    }
    if (transposed)
        vcf = false; // so that - is a transposed matrix on stdin
    if (vcf && alphabet_size == 0)
        alphabet_size = 2; // a site with more alleles stops the run, saying which -a to use
    if (minimal_block_size < 2)
    {
        std::cerr << "minimal_block_size is too small\n";
        exit(0);
    }
    if (alphabet_size < 2)
    {
        std::cerr << "alphabet_size is too small\n";
        exit(0);
    }

    std::cerr << "Running with alphabet: " << std::to_string(alphabet_size)
              << "\nwith minimal blocksize: " << std::to_string(minimal_block_size)
              << "\non file: " << filename << "\n";
    if (buffer_size)
        std::cerr << "with buffer size: " << std::to_string(buffer_size) << "\n";

    std::unique_ptr<ColumnReader> file;
    if (vcf)
    {
#ifdef WITH_HTSLIB
        if (buffer_size)
            std::cerr << "-g is ignored, a VCF is read one record at a time\n";
        file.reset(new VcfReader(filename, alphabet_size));
#else
        std::cerr << "This build cannot read VCF/BCF: build it with htslib (see the README)\n";
        exit(1);
#endif
    }
    else if (transposed)
    {
        if (buffer_size)
            std::cerr << "-g is ignored, a transposed matrix is read one column at a time\n";
        file.reset(new TransposedReader(filename, alphabet_size));
    }
    else
    {
        file.reset(new MatrixReader(filename, alphabet_size, buffer_size));
    }

    if (!file->next_column())
    {
        std::cerr << "No columns in the file: " << filename << "\n";
        exit(1);
    }
    PbwtOrder pbwt(file->column(), file->row_count(), alphabet_size);
    int_t columns = file->column_count(); // 0 if the file has to be read to know it
    int_t ak_max = 0, j = 0;
    for (bool more = true; more; j++)
    {
        if (j % (columns ? columns / 1000 + 1 : 1000) == 0)
        {
            std::cerr << "\rk:" << j << " #blocks: " << pbwt.total_blocks << " #SIZE: " << pbwt.get_current_size() << " #maxSize: " << ak_max;
            if (columns)
                std::cerr << " percent: " << (j * 100) / columns << "%";
            std::cerr << std::flush;
        }
        ak_max = std::max(ak_max, pbwt.get_current_size());

        pbwt.next();
        more = file->next_column();
        pbwt.gap_blocks(more);
    }
    std::cerr << "\n";
    std::cerr << "rows\t\t\t" << file->row_count() << "\n";
    std::cerr << "columns\t\t\t" << j << "\n";
    std::cerr << "ak_final_size\t\t" << pbwt.get_current_size() << "\n";
    std::cerr << "total_blocks_found\t" << pbwt.total_blocks << "\n";
    std::cerr << "total_expansions\t" << pbwt.expansion_count << "\n";
    std::cerr << "total_rows_collapsed\t" << pbwt.collapsed_rows_count << "\n";
    std::cerr << "total_range_collapsed\t" << pbwt.collapse_count << "\n";
}
