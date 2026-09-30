#ifndef __MATRIX_READER_HPP__
#define __MATRIX_READER_HPP__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef uint8_t allele_t;
const allele_t WILD = allele_t('*' - '0');

// class to read a matrix (one row per line) column by column
// buffer_size == 0: mmap, the whole file ends up in the page cache
// buffer_size > 0: one buffered stream per row, for panels that outgrow RAM (needs one open file per row)
class MatrixReader
{
  std::vector<allele_t> col;
  size_t rows;
  size_t cols;
  size_t size;
  char *mat = NULL;
  std::vector<std::ifstream> streams;
  std::vector<std::vector<char>> buffers;

public:
  MatrixReader(const std::string &filename, size_t buffer_size = 0)
  {
    auto fail = [&](const char *why)
    {
      std::cerr << why << ": " << filename << "\n";
      exit(1);
    };
    struct stat st;
    std::ifstream ifs(filename);
    std::string s;
    if (stat(filename.c_str(), &st) || !std::getline(ifs, s))
      fail("Couldn't open file");
    size = st.st_size;
    cols = s.size();
    rows = size / (cols + 1) + ((size % (cols + 1) != 0) ? 1 : 0);
    assert(size % (cols + 1) == 0 ||
           size % (cols + 1) == cols); // might be missing last eol
    col.resize(rows);

    if (buffer_size == 0)
    {
      int fd = open(filename.c_str(), O_RDONLY);
      mat = static_cast<char *>(mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0));
      close(fd);
      if (mat == MAP_FAILED)
        fail("Couldn't mmap file");
      return;
    }

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
  size_t getColSize() const { return cols; }
  size_t getRowSize() const { return rows; }

  // with streams the columns come in file order, whatever c is
  const allele_t *getColumn(size_t c)
  {
    for (size_t i = 0; i < rows; i++)
      col[i] = (mat ? mat[c + i * (cols + 1)] : streams[i].rdbuf()->sbumpc()) - '0';
    return col.data();
  }
};
#endif //__MATRIX_READER_HPP__
