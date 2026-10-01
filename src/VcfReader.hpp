#ifndef __VCF_READER_HPP__
#define __VCF_READER_HPP__
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "ColumnReader.hpp"

// The genotypes (FORMAT/GT) of a VCF or BCF, plain or compressed, or "-" for the standard input,
// one record (a site, so a column) at a time and in one pass. Haplotype h of sample s is the row
// s * ploidy + h, and a missing allele is a wildcard. Every record is a column, whatever its alleles:
// an allele above the alphabet size is an error (bcftools view -m2 -M2 -v snps keeps biallelic SNPs).
class VcfReader : public ColumnReader
{
  htsFile *fp;
  bcf_hdr_t *hdr;
  bcf1_t *rec;
  int32_t *gt = NULL;
  int gt_size = 0;
  size_t rows = 0;
  size_t records = 0;
  allele_t alphabet_size;
  bool unphased = false;
  std::vector<allele_t> col;

  void fail(const std::string &why)
  {
    std::cerr << why << " (record " << records << ")\n";
    exit(1);
  }

public:
  VcfReader(const std::string &filename, allele_t alphabet_size) : alphabet_size(alphabet_size)
  {
    fp = hts_open(filename.c_str(), "r");
    hdr = fp ? bcf_hdr_read(fp) : NULL;
    if (!hdr || bcf_hdr_nsamples(hdr) == 0)
    {
      std::cerr << "Couldn't read samples from file: " << filename << "\n";
      exit(1);
    }
    rec = bcf_init();
  }
  ~VcfReader()
  {
    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
  }
  size_t row_count() const { return rows; }
  size_t sample_count() const { return bcf_hdr_nsamples(hdr); }
  const allele_t *column() const { return col.data(); }

  bool next_column()
  {
    int r = bcf_read(fp, hdr, rec);
    if (r < -1)
      fail("Error reading the file");
    if (r == -1)
      return false;
    records++;

    int n = bcf_get_genotypes(hdr, rec, &gt, &gt_size);
    if (n <= 0)
      fail("No genotypes (FORMAT/GT)");
    if (rows == 0)
    {
      rows = n;
      col.resize(rows);
    }
    else if (size_t(n) != rows)
      fail("The number of haplotypes is not the one of the first record (ploidy changed?)");
    int ploidy = n / bcf_hdr_nsamples(hdr);

    for (int i = 0; i < n; i++)
    {
      int32_t g = gt[i];
      if (g == bcf_int32_vector_end || bcf_gt_is_missing(g)) // a lower ploidy is filled with the end marker
      {
        col[i] = WILD;
        continue;
      }
      int a = bcf_gt_allele(g);
      if (a >= alphabet_size)
        fail("Allele " + std::to_string(a) + " with alphabet size " + std::to_string(alphabet_size) + ", use -a " + std::to_string(a + 1));
      col[i] = a;
      if (i % ploidy && !bcf_gt_is_phased(g) && !unphased && a != bcf_gt_allele(gt[i - 1]))
      {
        unphased = true;
        std::cerr << "Warning: unphased genotypes (first at record " << records << "), taken in the order they are written\n";
      }
    }
    return true;
  }
};
#endif //__VCF_READER_HPP__
