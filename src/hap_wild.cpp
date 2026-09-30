#include <fstream>
#include <iostream>
#include <random>

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " <input_matrix> <wild_rate> <path_to_output_matrix>" << std::endl;
        std::cout << "Example: " << argv[0] << " data/input_matrix 3 data/output_matrix" << std::endl;

        return 1;
    }

    std::string filename = argv[1];
    double error = std::stod(argv[2]) / 100.0;
    std::string filename_errors = argv[3];

    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Couldn't open file: " << filename << "\n";
        return 1;
    }
    std::ofstream ofs(filename_errors);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> err(0, 1);

    long long int count = 0;
    std::string row;
    while (ifs >> row)
    {
        for (char &c : row)
        {
            if (err(gen) < error)
            {
                c = '*';
                count++;
            }
        }
        ofs << row << "\n";
    }

    std::cout << "added gaps: " << count << std::endl;
    std::cout << "filename: " << std::endl
              << filename << std::endl
              << "filename_errors: " << std::endl
              << filename_errors << std::endl;
    return 0;
}
