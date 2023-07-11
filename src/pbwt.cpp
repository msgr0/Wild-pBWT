#include "FileReader.hpp"
#include <iostream>
#include <vector>
#include <assert.h>
#include <sdsl/bit_vectors.hpp>
#include <unordered_set>
#include <ostream>
using namespace sdsl;

typedef uint8_t allele_t;
typedef long long int int_t;

// globals
bool verbose = false;
bool count_blocks = false;
bool output_blocks = false;
int_t minimal_block_size = 2;

void see(const std::vector<int_t> m)
{
    for (auto &i : m)
    {
        std::cerr << (int(i)) << " ";
    }
    std::cerr << "\n";
}

class PbwtOrder
{
private:
    std::vector<int_t> ak;
    std::vector<int_t> dk;
    std::vector<int_t> dk_sup;
    std::vector<bit_vector> b;
    bit_vector bv;

    int_t k;
    int_t N;
    int_t M;
    allele_t *column_pointer;

    allele_t alphabet_size;
    allele_t allele;
    int_t total_blocks;

    int_t expansion_count;
    int_t collapsed_rows_count;
    int_t collapse_count;
    int_t removed_count;

public:
    PbwtOrder(allele_t *r, int_t lines, int_t columns, allele_t alphabet_size)
    {
        k = 0;
        this->alphabet_size = alphabet_size;
        N = columns;
        M = lines;
        column_pointer = r;

        this->ak = std::vector<int_t>(this->M);
        for (int_t i = 0; i < this->M; i++)
        {
            ak[i] = i;
        }
        this->dk = std::vector<int_t>(this->M, 0);

        this->total_blocks = 0;
        this->expansion_count = 0;
        this->collapsed_rows_count = 0;
        this->collapse_count = 0;
    }
    void see()
    {
        std::cerr << "k: " << k << "\n";
        std::cerr << "column: \t";
        for (int_t i = 0; i < M; i++)
        {
            // std::cout << i << "\n";
            if (int(column_pointer[i]) == 250)
            {
                std::cerr << '*' << " ";
            }
            else
            {
                std::cerr << int(column_pointer[i]) << " ";
            }
        }
        std::cerr << "\n";

        std::cerr << "ordered: \t";
        for (auto &i : ak)
        {
            if (int(column_pointer[i]) == 250)
            {
                std::cerr << '*' << " ";
            }
            else
            {
                std::cerr << int(column_pointer[i]) << " ";
            }
        }
        std::cerr << "\n";

        std::cerr << "ak: \t\t";
        for (auto &i : ak)
        {
            std::cerr << int(i) << " ";
        }
        std::cerr << "\n";

        std::cerr << "dk: \t\t";
        for (auto &i : dk)
        {
            std::cerr << int(i) << " ";
        }
        std::cerr << "\n";

        std::cerr << "bitvectors: \t";
        for (auto &i : b)
        {
            std::cerr << "b[" << i << "] ";
            std::cerr << "\n";
        }
        std::cerr << "\n";
        std::cerr << "------------------------\n";
    }

    int_t get_total_blocks()
    {
        return total_blocks;
    }
    int_t get_expansion_count()
    {
        return expansion_count;
    }
    int_t get_collapsed_rows_count()
    {
        return collapsed_rows_count;
    }
    int_t get_collapse_count()
    {
        return collapse_count;
    }

    std::vector<int_t> get_curr_ak()
    {
        return ak;
    }

    std::vector<int_t> get_curr_dk() { return dk; }

    int_t get_current_k() { return k; }
    int_t get_current_size() { return ak.size(); }

    int_t get_N() { return N; }

    void gap_blocks()
    {
        int_t current_size = get_current_size();
        allele_t next_allele = 0;

        this->b = std::vector<bit_vector>(alphabet_size);
        for (int_t i = 0; i < alphabet_size; i++)
        {
            b[i] = bit_vector(current_size, 0);
        }
        this->dk_sup = std::vector<int_t>(current_size, -1);
        for (int_t i = 0; i < current_size; i++)
        {
            if (k >= N - 1)
            {
                for (int_t l = 0; l < alphabet_size; l++)
                {
                    b[l][i] = 0;
                }
            }
            else
            {
                next_allele = column_pointer[ak[i]];
                if (next_allele == 250)
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
        std::vector<rank_support_v<0>> rank(alphabet_size + 1);
        for (int_t l; l < alphabet_size; l++)
        {
            rank[l] = rank_support_v<0>(&b[l]);
        }
        if (k >= 0 && k < N)
        {
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
                        if (next_allele == 250)
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
                        if (next_allele == 250)
                        {
                            next_allele = column_pointer[ak[n]];
                        }
                        if (ak[n] != ak[n - 1])
                            skip = false;
                        n++;
                    }
                    n--;

                    if (next_allele != 250 && !skip)
                    {
                        int_t width = k - dk[i] + 1;
                        int_t diff;
                        diff = rank[next_allele](n + 1) - rank[next_allele](m);
                        if (diff != 0 || k >= N - 1)
                        {
                            dk_sup[m] = dk[i];
                            if (count_blocks)
                            {
                                total_blocks += 1;
                            }
                            else
                            {
                                std::unordered_set<int_t> b;
                                for (int_t a = m; a <= n; a++)
                                {
                                    b.insert(ak[a]);
                                }
                                int_t range = b.size();

                                if ((width)*range >= minimal_block_size)
                                {
                                    total_blocks++;
                                    if (output_blocks)
                                    {
//
                                    //    std::cout
                                    //        << range * width << " (" << dk[i] << "-"
                                    //        << k << "," << ak[m] << ":" << range << ")\n"; // \t{";
std::cout << "[";
                                        for (auto &b : b)
                                        {
                                            std::cout << b << ",";
                                        }
std::cout << "], " << dk[i] << ", " << k << "\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        this->k++;
    }

    void next(allele_t *r)
    {
        int_t size_curr = 0;
        std::vector<std::vector<int_t>> a(alphabet_size, std::vector<int_t>(0));
        std::vector<std::vector<int_t>> d(alphabet_size, std::vector<int_t>(0));
        std::vector<int_t> p(alphabet_size, 0);
        std::vector<int_t> u(alphabet_size, 0);
        std::vector<std::vector<int_t>> m(alphabet_size, std::vector<int_t>(0));

        for (auto &i : a)
        {
            i.reserve(get_current_size());
        }

        for (auto &i : d)
        {
            i.reserve(get_current_size());
        }

        for (int_t l = 0; l < alphabet_size; l++)
        {
            u[l] = 0;
            p[l] = k + 1;
        }

        for (int_t i = 0; i < get_current_size(); i++)
        {
            allele = column_pointer[ak[i]];

            for (int_t l = 0; l < alphabet_size; l++)
            {
                if (dk[i] > p[l])
                {
                    p[l] = dk[i];
                }
            }
            if (allele != 250)
            {
                a[allele].emplace_back(ak[i]);
                d[allele].emplace_back(p[allele]);
                p[allele] = 0;
                u[allele] = u[allele] + 1;
                size_curr += 1;
            }
            else
            {
                // assert(allele == 255);
                this->expansion_count += 1;
                for (int_t m = 0; m < alphabet_size; m++)
                {
                    a[m].emplace_back(ak[i]);
                    d[m].emplace_back(p[m]);
                    p[m] = 0;
                    u[m] = u[m] + 1;
                    size_curr += 1;
                }
            }
        }

        ak.clear();
        dk.clear();
        ak.reserve(size_curr);
        dk.reserve(size_curr);

        for (int_t i = 0; i < alphabet_size; i++)
        {
            for (int_t j = 0; j < a[i].size(); j++)
                ak.emplace_back(a[i][j]);
            for (int_t j = 0; j < d[i].size(); j++)
                dk.emplace_back(d[i][j]);
        }

        std::vector<int_t> ac(0);
        std::vector<int_t> dc(0);
        ac.reserve(size_curr);
        dc.reserve(size_curr);

        bool started = false;
        int_t max_d = 0;
        int_t upper_d = 0;
        int_t begin = 0;
        for (int_t i = 0; i < size_curr - 1; i++)
        {
            if (ak[i] == ak[i + 1])
            {
                if (!started)
                {
                    started = true;
                    upper_d = dk[i];
                    begin = ak[i];
                    ac.emplace_back(begin);
                    dc.emplace_back(upper_d);
                    max_d = max(dk[i + 1], max_d);
                }
                else
                {
                    max_d = max(dk[i + 1], max_d);
                    this->collapsed_rows_count += 1;
                    removed_count += 1;
                }
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

        ak.clear();
        dk.clear();
        this->ak = ac;
        this->dk = dc;

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
    std::cerr << "-b <min block size> -g <buffer_size:512> -v y <verbose> \n";
    std::cerr << "Example: ./bin/wild-pbwt -a 3 -f paper_wild -o y > paper_wild_out.txt \n";
    exit(0);
}

int main(int argc, char **argv)
{
    int ch;
    std::string filename;
    allele_t alphabet_size;
    int_t buffer_size = 512;

    while ((ch = getopt(argc, argv, "h:f:a:c:o:v:g:b:")) != -1)
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
    if (buffer_size < 2)
    {
        std::cerr << "buffersize is too small\n";
        exit(0);
    }
    if (alphabet_size < 2)
    {
        std::cerr << "alphabet_size is too small\n";
        exit(0);
    }

    // bool gaps = (with_gaps == 0) ? false : true;
    std::cerr << "Running with alphabet: " << std::to_string(alphabet_size)
              << "\nwith minimal blocksize: " << std::to_string(minimal_block_size)
              << "\non file: " << filename
              << "\nwith buffer size: " << std::to_string(buffer_size) << "\n";

    LR_file_hap file(filename, buffer_size);

    int_t columns = file.get_number_column();
    int_t lines = file.get_number_line();
    allele_t *ni = file.next();
    int_t j = 0;
    PbwtOrder pbwt = PbwtOrder(ni, lines, columns, alphabet_size);
    size_t ak_max = 0;
    while (!file.is_end())
    {

        if (j % (int(columns / 1000) + 1) == 0)
        {
            std::cerr << "\rk:" << j << " #blocks: " << pbwt.get_total_blocks() << " #SIZE: " << pbwt.get_curr_ak().size() << " #maxSize: " << ak_max << " percent: " << (j * 100) / columns << "%" << std::flush;
        }
        ak_max = std::max(ak_max, pbwt.get_curr_ak().size());
        j++;

        pbwt.next(ni);
        ni = file.next();
        pbwt.gap_blocks();
    }
    std::cerr << "\n";
    std::cerr << "ak_final_size\t\t" << pbwt.get_curr_ak().size() << "\n";
    std::cerr << "total_blocks_found\t" << pbwt.get_total_blocks() << "\n";
    std::cerr << "total_expansions\t" << pbwt.get_expansion_count() << "\n";
    std::cerr << "total_rows_collapsed\t" << pbwt.get_collapsed_rows_count() << "\n";
    std::cerr << "total_range_collapsed\t" << pbwt.get_collapse_count() << "\n";
}
