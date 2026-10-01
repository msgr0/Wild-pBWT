#ifndef __MATRIX_READER_HPP__
#define __MATRIX_READER_HPP__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "ColumnReader.hpp"

// class to read a matrix (one row per line, as many characters in each) column by column
// buffer_size == 0: mmap, the whole file ends up in the page cache
// buffer_size > 0: one buffered stream per row, for panels that outgrow RAM (needs one open file per row)
class MatrixReader : public ColumnReader
{
  std::vector<allele_t> col;
  size_t rows;
  size_t cols;
  size_t size;
  size_t next = 0;
  allele_t alphabet_size;
  std::string name;
  char *mat = NULL;
  std::vector<std::ifstream> streams;
  std::vector<std::vector<char>> buffers;

  void fail(const std::string &why)
  {
    std::cerr << why << ": " << name << "\n";
    exit(1);
  }

public:
  MatrixReader(const std::string &filename, allele_t alphabet_size, size_t buffer_size = 0)
      : alphabet_size(alphabet_size), name(filename)
  {
    struct stat st;
    std::ifstream ifs(filename);
    std::string s;
    if (stat(filename.c_str(), &st) || !ifs)
      fail("Couldn't open file");
    if (!std::getline(ifs, s) || s.empty())
      fail("The file, or its first line, is empty");
    if (s.back() == '\r')
      fail("The lines end with \\r (a Windows file?), remove it for example with tr -d '\\r'");
    size = st.st_size;
    cols = s.size();
    rows = (size + cols) / (cols + 1); // the lines with a newline, and one without
    col.resize(rows);

    if (buffer_size == 0)
    {
      int fd = open(filename.c_str(), O_RDONLY);
      mat = static_cast<char *>(mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0));
      close(fd);
      if (mat == MAP_FAILED)
        fail("Couldn't mmap file");
    }

    // The columns are found by the position, so every line must have as many characters as the first
    // (and the reads stay inside the file): line i has to end where it would if all were the same.
    for (size_t i = 0; i < rows; i++)
    {
      size_t end = (cols + 1) * i + cols;
      int c = EOF;
      if (end < size)
      {
        if (mat)
          c = mat[end];
        else
          c = (ifs.clear(), ifs.seekg(end), ifs.get());
      }
      if (end > size || !(c == '\n' || (i == rows - 1 && end == size)))
        fail("Line " + std::to_string(i + 1) + " does not have " + std::to_string(cols) + " characters, as the first one has");
    }
    if (mat)
      return;

    struct rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    limit.rlim_cur = std::min<rlim_t>(limit.rlim_max, std::max<rlim_t>(limit.rlim_cur, rows + 64));
    setrlimit(RLIMIT_NOFILE, &limit);

    streams = std::vector<std::ifstream>(rows);
    buffers.assign(rows, std::vector<char>(buffer_size));
    for (size_t i = 0; i < rows; i++)
    {
      streams[i].rdbuf()->pubsetbuf(buffers[i].data(), buffer_size);
      streams[i].open(filename);
      streams[i].seekg((cols + 1) * i);
      if (!streams[i])
        fail("Couldn't open one stream per row (open files limit?), retry without -g");
    }
  }
  ~MatrixReader()
  {
    if (mat)
      munmap(mat, size);
  }
  size_t row_count() const { return rows; }
  size_t column_count() const { return cols; }
  const allele_t *column() const { return col.data(); }

  bool next_column()
  {
    if (next >= cols)
      return false;
    for (size_t i = 0; i < rows; i++)
    {
      char c = mat ? mat[next + i * (cols + 1)] : streams[i].rdbuf()->sbumpc();
      col[i] = read_allele(c, alphabet_size);
      if (col[i] == NOT_ALLELE)
        fail(bad_allele(c, alphabet_size) + " in haplotype " + std::to_string(i) + ", column " + std::to_string(next));
    }
    next++;
    return true;
  }
};
#endif //__MATRIX_READER_HPP__
