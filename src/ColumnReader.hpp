#ifndef __COLUMN_READER_HPP__
#define __COLUMN_READER_HPP__
#include <cstddef>
#include <cstdint>
#include <string>

typedef uint8_t allele_t;
const allele_t WILD = allele_t('*' - '0');
const allele_t NOT_ALLELE = 255;

// The allele that a text matrix writes as the character c: the digit 0, 1, ... below the alphabet
// size, or WILD for *. Anything else is NOT_ALLELE.
inline allele_t read_allele(char c, allele_t alphabet_size)
{
  int a = c - '0';
  if (c == '*')
    return WILD;
  return a >= 0 && a < alphabet_size ? allele_t(a) : NOT_ALLELE;
}

// What is wrong with a character that read_allele refused.
inline std::string bad_allele(char c, allele_t alphabet_size)
{
  int a = c - '0';
  if (a >= alphabet_size && a < 10)
    return "Allele " + std::to_string(a) + " with alphabet size " + std::to_string(alphabet_size) + ", use -a " + std::to_string(a + 1);
  unsigned char u = c;
  if (u >= 32 && u < 127)
    return std::string("Unexpected character '") + c + "'";
  const char *hex = "0123456789abcdef";
  return std::string("Unexpected character 0x") + hex[u >> 4] + hex[u & 15];
}

// A panel read one column (one site) at a time. next_column() puts the alleles of the next site
// in column(), one per haplotype and WILD for a wildcard, and tells whether there was one.
// The pointer returned by column() stays the same, and row_count() is known after the first column.
class ColumnReader
{
public:
  virtual ~ColumnReader() {}
  virtual size_t row_count() const = 0;
  virtual size_t column_count() const { return 0; } // 0 if not known before reading everything
  virtual const allele_t *column() const = 0;
  virtual bool next_column() = 0;
};
#endif
