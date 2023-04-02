#ifndef __MATRIX_READER_HPP__
#define __MATRIX_READER_HPP__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>

typedef std::vector<uint8_t> matrixType;

class MatrixReader
{
public:
  enum Method
  {
    M_byCol, // assumes matrix is given columnwise
    M_mmap   // read column on the fly by mmap
  };
  size_t currentCol;
  bool byCol;
  std::ifstream ifs;
  std::vector<matrixType> matrix;
  matrixType nextCol;
  size_t rows;
  size_t cols;
  enum Method method;

  struct stat st;
  int fd;
  char *mat;

public:
  // class to read matrix from file
  // givenByCol = true: assumes that the file is transposed
  MatrixReader(const std::string &filename, enum Method m)
      : currentCol(0), method(m)
  {
    if (stat(filename.c_str(), &st))
    {
      std::cerr << "Couldn't open file: " << filename << "\n";
      exit(1);
    }
    std::string s;
    ifs.open(filename, std::ifstream::in);
    ifs >> s;
    ifs.close();
    cols = s.size();
    rows = st.st_size / (cols + 1) + ((st.st_size % (cols + 1) != 0) ? 1 : 0);
    assert(st.st_size % (cols + 1) == 0 ||
           st.st_size % (cols + 1) == cols); // might be missing last eol

    switch (method)
    {
    case M_byCol:
      ifs.open(filename, std::ifstream::in);
      break;
    case M_mmap:
      fd = open(filename.c_str(), O_RDONLY);
      mat = static_cast<char *>(
          mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0));
      nextCol.resize(rows);
      break;
    }
  };
  size_t getColSize() const { return cols; }
  size_t getRowSize() const { return rows; }

  const matrixType &getNextColumn()
  {
    std::string s;
    switch (method)
    {
    case M_byCol:
      ifs >> s;
      nextCol.resize(s.size());
      for (size_t i = 0; i < s.size(); i++)
        nextCol[i] = (s[i] == '1');
      if (s.empty())
        ifs.close();
      return nextCol;
      break;

    case M_mmap:
      if (currentCol < cols)
      {
        size_t off = currentCol;
        for (size_t i = 0; i < rows; i++)
        {
          if (mat[off] == '*')
          {
            nextCol[i] = -1; // 255 uint
          }
          else
          {
            nextCol[i] = mat[off] - 48;
          }
          off += cols + 1;
        }
        currentCol++;
        return nextCol;
      }
      else
      {
        nextCol.clear();
        return nextCol;
      }
      return nextCol;
      break;
    }
    return nextCol; // WARNING ==>  control reaches end of non-void function
  };
  ~MatrixReader()
  {
    switch (method)
    {
    case M_mmap:
      munmap(static_cast<void *>(mat), st.st_size);
      close(fd);
      break;
    default:
      break;
    }
  }
};
#endif //__MATRIX_READER_HPP__
