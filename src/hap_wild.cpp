#include "./MatrixReader.hpp"
#include <iostream>
#include <ostream>
#include <random>
#include <iostream>
#include <fstream>
#include <thread>

int main(int argc, char **argv)
{
    int alphabet = 2;
    std::fstream ofs;

    long long int count = 0;

    int seed = 0;
    seed = seed == 0 ? std::random_device()() : seed;
    std::mt19937 gen(seed);

    if (argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " <input_matrix> <wild_rate> <path_to_output_matrix>" << std::endl;
        std::cout << "Example: " << argv[0] << " data/input_matrix 3 data/output_matrix" << std::endl;

        return 1;
    }

    std::string save_dir = argv[3];

    double error = std::stod(argv[2]) / 100.0;
    std::string filename = argv[1];

    std::uniform_real_distribution<> err(0, 1);

    MatrixReader::Method method = MatrixReader::M_byCol;
    MatrixReader matrix(filename, method);
    int n = matrix.getColSize();
    int m = matrix.getRowSize();

    ofs.open(save_dir, std::ofstream::out);

    for (int i = 0; i < m; i++)
    {
        auto column = matrix.getNextColumn();
        for (int j = 0; j < n; j++)
        {
            if (err(gen) < error)
            {
                ofs << "*";
                count++;
            }
            else
                ofs << column[j] + 0;
        }

        ofs << std::endl;
        std::cout << "\rpourcent: " << (i * 100) / m << "%" << std::flush;
    }
    ofs.close();

    std::cout << std::endl
              << "added gaps: " << count << std::endl;
    std::cout << "filename: " << std::endl
              << filename << std::endl
              << "filename_errors: " << std::endl
              << save_dir << std::endl;
    return 0;
}
