#ifndef __FILE_READER_HPP__
#define __FILE_READER_HPP__

#define LIMITSIZE 20000

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <cstdint>

#include <iomanip>

// #include <iostream>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <algorithm>

#include <tuple>
#include <list>

#include <unistd.h>

using namespace std;

typedef uint8_t allele_t;
typedef long long int int_t;

class LR_file_hap
{
private:
    // Attributs
    string filename;

    string alphabet;
    int_t *alphabet_int;
    int_t alphabet_size;

    int_t nb_line;
    int_t nb_column;
    int_t file_size;
    int_t current_column;
    ifstream *filearray;
    istreambuf_iterator<char> *itarray;
    char **filebuffer;

    allele_t *current_char;

public:
    int_t get_number_line(void)
    {
        return nb_line;
    }

    int_t get_number_column(void)
    {
        return nb_column;
    }

    bool is_end(void)
    {
        return current_column > nb_column;
    }

    void close(void)
    {
        for (int_t i = 0; i < get_number_line(); i++)
        {
            filearray[i].close();
        }
    }

    allele_t *next(void)
    {
        if (current_column > nb_column)
        {
            return NULL;
        }
        else
        {
            for (int_t i = 0; i < nb_line; i++)
            {
                current_char[i] = *itarray[i]++ - '0';
            }
            current_column++;
            return current_char;
        }
    }

    allele_t *get_current_char(void)
    {
        return current_char;
    }

    LR_file_hap(string filename, int_t buffer_size)
    {

        struct rlimit limit;

        limit.rlim_cur = LIMITSIZE;
        limit.rlim_max = LIMITSIZE;

        if (setrlimit(RLIMIT_NOFILE, &limit) != 0)
        {
            cerr << "setrlimit() failed with errno=" << errno << endl;
        }

        if (getrlimit(RLIMIT_NOFILE, &limit) != 0)
        {
            cerr << "getrlimit() failed with errno=" << errno << endl;
        }

        current_column = 0;
        nb_line = 0;
        nb_column = 0;
        file_size = 0;

        // int BufferSize = 1024;
        int BufferSize_big = 8192;
        int BufferSize = buffer_size;

        ifstream str2;
        char _buffer[BufferSize_big];
        str2.rdbuf()->pubsetbuf(_buffer, BufferSize_big);
        str2.open(filename.c_str(), ifstream::in);
        if (!str2.is_open())
        {
            cerr << "Error: file " << filename << " not found" << endl;
            exit(1);
        }
        istreambuf_iterator<char> istart(str2);
        istreambuf_iterator<char> iend;

        while (istart != iend)
        {
            if (*istart == '\n')
            {
                nb_line++;
            }
            file_size++;
            istart++;
        }

        nb_column = (int_t)(file_size / nb_line - 1);

        filearray = new ifstream[nb_line];
        filebuffer = new char *[nb_line];
        itarray = new istreambuf_iterator<char>[nb_line];
        //
        str2.clear();
        str2.seekg(0);
        istart = str2;
        //
        //
        istreambuf_iterator<char> it;

        for (int_t i = 0; i < nb_line; i++)
        {
            filebuffer[i] = new char[BufferSize];
            filearray[i].rdbuf()->pubsetbuf(filebuffer[i], BufferSize);

            filearray[i].open(filename.c_str(), ifstream::in);

            filearray[i].seekg((nb_column + 1) * i, ios::beg);
            itarray[i] = filearray[i];
        }

        current_char = new allele_t[nb_line];

        str2.close();
    }
};

#endif //__FILE_READER_HPP__
