#ifndef __TRANSPOSED_READER_HPP__
#define __TRANSPOSED_READER_HPP__
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "ColumnReader.hpp"

// A matrix written by columns: every line is a site with one character per haplotype (0, 1, ... and *
// for a wildcard), so the file is the transpose of the usual one, and it is read in one pass (also
// from the standard input, "-"). Empty lines are skipped, a \r at the end of a line is ignored.
class TransposedReader : public ColumnReader
{
  std::ifstream file;
  std::istream *in;
  std::string line;
  std::vector<allele_t> col;
  size_t rows = 0;
  size_t line_number = 0;
  allele_t alphabet_size;

  void fail(const std::string &why)
  {
    std::cerr << why << " (line " << line_number << ")\n";
    exit(1);
  }

public:
  TransposedReader(const std::string &filename, allele_t alphabet_size) : alphabet_size(alphabet_size)
  {
    if (filename == "-")
    {
      in = &std::cin;
      return;
    }
    file.open(filename);
    if (!file)
    {
      std::cerr << "Couldn't open file: " << filename << "\n";
      exit(1);
    }
    in = &file;
  }
  size_t row_count() const { return rows; }
  const allele_t *column() const { return col.data(); }

  bool next_column()
  {
    do
    {
      if (!std::getline(*in, line))
        return false;
      line_number++;
      if (!line.empty() && line.back() == '\r')
        line.pop_back();
    } while (line.empty());

    if (rows == 0)
    {
      rows = line.size();
      col.resize(rows);
    }
    else if (line.size() != rows)
      fail("The line has " + std::to_string(line.size()) + " characters, the first one has " + std::to_string(rows));

    for (size_t i = 0; i < rows; i++)
    {
      col[i] = read_allele(line[i], alphabet_size);
      if (col[i] == NOT_ALLELE)
        fail(bad_allele(line[i], alphabet_size) + " in haplotype " + std::to_string(i));
    }
    return true;
  }
};
#endif //__TRANSPOSED_READER_HPP__
