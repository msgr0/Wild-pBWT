#include <fstream>
#include <iostream>
#include <random>

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        std::cout << "Usage: " << argv[0] << " <save_directory> <t-alleles> <haplotypes_count> <SNPs_count>" << std::endl;
        std::cout << "Example: " << argv[0] << " data 3 5000 100000" << std::endl;

        return 1;
    }
    std::string save_dir = argv[1];
    save_dir.append("/");
    int alphabet = std::stoi(argv[2]);
    int m = std::stoi(argv[3]);
    int n = std::stoi(argv[4]);

    std::string filename = "hap" + std::to_string(alphabet) + "_gen_" + std::to_string(m) + "_" + std::to_string(n);
    std::ofstream ofs(save_dir + filename);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis(0, alphabet - 1);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ofs << dis(gen);
        }

        ofs << std::endl;
        std::cout << "\rpourcent: " << (i * 100) / m << "%" << std::flush;
    }
    std::cout << std::endl
              << "filename: " << std::endl
              << filename << std::endl;
    return 0;
}
