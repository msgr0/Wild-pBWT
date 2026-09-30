#include "MatrixReader.hpp"
#include <iostream>
#include <numeric>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <unordered_set>
using namespace sdsl;

typedef long long int int_t;

// globals
bool verbose = false;
bool count_blocks = false;
bool output_blocks = false;
int_t minimal_block_size = 2;

class PbwtOrder
{
private:
    std::vector<int_t> ak;
    std::vector<int_t> dk;
    std::vector<int_t> dk_sup;
    std::vector<bit_vector> b;
    std::vector<int_t> in_block; // per row, last block/group it was seen in (see maximal_rows)
    std::vector<int_t> in_group;
    int_t stamp = 0;
    std::vector<std::vector<int_t>> where; // per row, the positions of its copies in ak

    int_t k = 0;
    int_t N;
    int_t M;
    const allele_t *column_pointer;

    allele_t alphabet_size;

public:
    int_t total_blocks = 0;
    int_t expansion_count = 0;
    int_t collapsed_rows_count = 0;
    int_t collapse_count = 0;

    PbwtOrder(const allele_t *r, int_t lines, int_t columns, allele_t alphabet_size)
        : ak(lines), dk(lines, 0), in_block(lines, 0), in_group(lines, 0), where(lines), N(columns), M(lines), column_pointer(r), alphabet_size(alphabet_size)
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
        std::cerr << "\nbitvectors: \t";
        for (auto &i : b)
            std::cerr << "b[" << i << "] \n";
        std::cerr << "\n------------------------\n";
    }

    int_t get_current_size() { return ak.size(); }

    // number of distinct rows in ak[m..n], or 0 if the block (rows, d, k) must not be
    // reported: it is not left-maximal, not row-maximal, or it is reported from elsewhere.
    //
    // Left: dk only says where two wildcard-filled copies start to match, so a divergence
    // at column d - 1 may come from two fillings of a wildcard and not from two rows that
    // really differ. The entries split, where dk == d, into groups whose copies agree on
    // column d - 1: if one group has all the rows, the rows are compatible on that column
    // and the block extends to the left.
    //
    // Rows: when all the rows have a wildcard on the same column, every other filling of it
    // gives another run of entries matching from d, with the same rows plus the ones that
    // really have that allele. Any such run holds a copy of each row of the block, so the
    // copies of one row find them all: a run with more rows means the block is not
    // row-maximal, a run with the same rows is the same block, reported from its first run.
    int_t maximal_rows(int_t m, int_t n, int_t d)
    {
        int_t block = ++stamp, group = ++stamp;
        int_t range = 0, group_range = 0, largest_group = 0;
        int_t fewest = ak[m]; // row of the block with the fewest copies
        for (int_t a = m; a <= n; a++)
        {
            if (a > m && dk[a] == d)
            {
                group = ++stamp;
                group_range = 0;
            }
            if (in_block[ak[a]] != block)
            {
                in_block[ak[a]] = block;
                range++;
                if (where[ak[a]].size() < where[fewest].size())
                    fewest = ak[a];
            }
            if (in_group[ak[a]] != group)
            {
                in_group[ak[a]] = group;
                largest_group = std::max(largest_group, ++group_range);
            }
        }
        if (largest_group == range)
            return 0;

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

    void gap_blocks()
    {
        int_t current_size = get_current_size();
        allele_t next_allele = 0;

        b.assign(alphabet_size, bit_vector(current_size, 0));
        dk_sup.assign(current_size, -1);
        for (auto &w : where)
        {
            w.clear();
        }
        for (int_t i = 0; i < current_size; i++)
        {
            where[ak[i]].push_back(i);
        }
        if (k < N - 1)
        {
            for (int_t i = 0; i < current_size; i++)
            {
                next_allele = column_pointer[ak[i]];
                if (next_allele == WILD)
                {
                    for (int_t l = 0; l < alphabet_size; l++)
                    {
                        b[l][i] = 1;
                    }
                }
                else
                {
                    b[next_allele][i] = 1;
                }
            }
        }
        std::vector<rank_support_v<0>> rank(alphabet_size);
        for (int_t l = 0; l < alphabet_size; l++)
        {
            rank[l] = rank_support_v<0>(&b[l]);
        }
        for (int_t i = 1; i < current_size; i++)
        {
            bool skip = true;
            if (k == N - 1)
                next_allele = 0;
            else
            {
                next_allele = column_pointer[ak[i]];
            }
            if (dk[i] <= k)
            {
                int_t m = i;
                int_t n = i + 1;

                while (m > 0 && dk[m] <= dk[i])
                {
                    if (ak[m] != ak[m - 1])
                        skip = false;
                    m--;
                    if (next_allele == WILD)
                    {
                        next_allele = column_pointer[ak[m]];
                    }
                }

                if (dk_sup[m] == dk[i])
                {
                    continue;
                }

                while (n < current_size && dk[n] <= dk[i])
                {
                    if (next_allele == WILD)
                    {
                        next_allele = column_pointer[ak[n]];
                    }
                    if (ak[n] != ak[n - 1])
                        skip = false;
                    n++;
                }
                n--;

                if (next_allele != WILD && !skip)
                {
                    int_t width = k - dk[i] + 1;
                    int_t diff = rank[next_allele](n + 1) - rank[next_allele](m);
                    if (diff != 0 || k >= N - 1)
                    {
                        dk_sup[m] = dk[i];
                        int_t range = maximal_rows(m, n, dk[i]);
                        if (range == 0)
                        {
                            continue;
                        }
                        if (count_blocks)
                        {
                            total_blocks += 1;
                        }
                        else if ((width)*range >= minimal_block_size)
                        {
                            total_blocks++;
                            if (output_blocks)
                            {
                                std::unordered_set<int_t> rows(ak.begin() + m, ak.begin() + n + 1);
                                std::cout << "[";
                                for (auto &r : rows)
                                {
                                    std::cout << r << ",";
                                }
                                std::cout << "], " << dk[i] << ", " << k << "\n";
                            }
                        }
                    }
                }
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
                    if (upper_d == max_d)
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
    std::cerr << "-a <alphabet_size> \n-f <filename> \n";
    std::cerr << "-c y <count max blocks> \n-o y <out_blocks to std_out> \n";
    std::cerr << "-b <min block size> -v y <verbose> \n";
    std::cerr << "-g <buffer_size> <read with one buffered stream per row instead of mmap> \n";
    std::cerr << "Example: ./bin/wild-pbwt -a 3 -f paper_wild -o y > paper_wild_out.txt \n";
    exit(0);
}

int main(int argc, char **argv)
{
    int ch;
    std::string filename;
    allele_t alphabet_size = 0;
    int_t buffer_size = 0;

    while ((ch = getopt(argc, argv, "hf:a:c:o:v:g:b:")) != -1)
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

    MatrixReader file(filename, buffer_size);

    int_t columns = file.getColSize();
    int_t lines = file.getRowSize();
    PbwtOrder pbwt = PbwtOrder(file.getColumn(0), lines, columns, alphabet_size);
    int_t ak_max = 0;
    for (int_t j = 0; j < columns; j++)
    {

        if (j % (int(columns / 1000) + 1) == 0)
        {
            std::cerr << "\rk:" << j << " #blocks: " << pbwt.total_blocks << " #SIZE: " << pbwt.get_current_size() << " #maxSize: " << ak_max << " percent: " << (j * 100) / columns << "%" << std::flush;
        }
        ak_max = std::max(ak_max, pbwt.get_current_size());

        pbwt.next();
        if (j + 1 < columns)
            file.getColumn(j + 1);
        pbwt.gap_blocks();
    }
    std::cerr << "\n";
    std::cerr << "ak_final_size\t\t" << pbwt.get_current_size() << "\n";
    std::cerr << "total_blocks_found\t" << pbwt.total_blocks << "\n";
    std::cerr << "total_expansions\t" << pbwt.expansion_count << "\n";
    std::cerr << "total_rows_collapsed\t" << pbwt.collapsed_rows_count << "\n";
    std::cerr << "total_range_collapsed\t" << pbwt.collapse_count << "\n";
}
