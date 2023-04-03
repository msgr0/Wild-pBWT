#include <random>
#include <iostream>
#include <fstream>
#include <thread>

int main(int argc, char **argv)
{
    if (argc != 5 && argc != 6)
    {
        // std::cout << "Usage: " << argv[0] << " <save_dir> <m> <n> <error_percentage> " << std::endl;
        std::cout << "Usage: " << argv[0] << " <save_directory> <t-alleles> <haplotypes_count> <SNPs_count> <wild_rate> " << std::endl;
        std::cout << "Example: " << argv[0] << " data 3 5000 100000 3" << std::endl;

        return 1;
    }
    std::string save_dir = argv[1];
    save_dir.append("/");
    int alphabet = std::stoi(argv[2]);
    int m = std::stoi(argv[3]);
    int n = std::stoi(argv[4]);

    std::string filename = "hap" + std::to_string(alphabet) + "_gen_" + std::to_string(m) + "_" + std::to_string(n);
    double error = 0;
    std::string filename_errors;
    if (argc == 6)
    {
        std::cout << "adding gaps" << std::endl;
        error = std::stod(argv[5]) / 100;
        filename_errors = filename + "_wild_" + std::to_string(error * 100).substr(0, 5);
    }
    std::fstream ofs;
    std::fstream ifs;
    std::fstream ofs_err;

    long long int count = 0;
    ofs.open(save_dir + filename, std::ofstream::out);
    ifs.open(save_dir + filename, std::ofstream::in);

    if (argc == 6)
        ofs_err.open(save_dir + filename_errors, std::ofstream::out);
    int seed = 0;
    seed = seed == 0 ? std::random_device()() : seed;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, alphabet - 1);
    std::uniform_real_distribution<> err(0, 1);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            auto add = dis(gen);
            if (argc == 6)
            {
                if (err(gen) < error)
                {
                    ofs_err << "*";
                    count++;
                }
                else
                    ofs_err << add;
            }

            ofs << add;
        }

        ofs << std::endl;
        ofs_err << std::endl;
        std::cout << "\rpourcent: " << (i * 100) / m << "%" << std::flush;
    }
    ofs.close();
    if (argc == 6)
        ofs_err.close();
    std::cout << std::endl
              << "added gaps: " << count << std::endl;
    std::cout << "filename: " << std::endl
              << filename << std::endl
              << "filename_errors: " << std::endl
              << filename_errors << std::endl;
    return 0;
}
