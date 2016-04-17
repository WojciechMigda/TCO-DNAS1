/*******************************************************************************
 * Copyright (c) 2016 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: dna_sequencing.hpp
 *
 * Description:
 *      description
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2016-04-16   wm              Initial version
 *
 ******************************************************************************/


#ifndef DNA_SEQUENCING_HPP_
#define DNA_SEQUENCING_HPP_

#include "BurrowsWheeler.hpp"
#include "gx/dna_nucleobase.hpp"
#include "gx/dna_nucleobase_istream.hpp"

#include "boost/algorithm/string/split.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>


struct DNASequencing
{
    void initTest()
    {
        m_curr_chromaid = -1;
        m_subchroma.clear();
    }

    int initTest(int testDifficulty)
    {
        initTest();
        return 0;
    }

    int preProcessing()
    {
        // close any sub chroma remaining after passReferenceGenome
        close_sub_chroma();

        return 0;
    }

    void close_sub_chroma()
    {
        constexpr std::size_t MIN_SUBCHROMA_SZ = 300;

        // process subchroma collected so far
        if (m_subchroma.size() >= MIN_SUBCHROMA_SZ)
        {
//            std::stringstream ss(m_subchroma);
//            std::vector<gx::DnaNucleobase> subchroma;
//
//            ss >> subchroma;
            // PY: first_occurence = FirstOccurence(Text)
            // FirstOccurence(text) -> indeksy pierwszych wystapien w posortowanym text
            // text = ACATGCTACTTT$
            // ordered = [u'$', u'A', u'A', u'A', u'C', u'C', u'C', u'G', u'T', u'T', u'T', u'T', u'T']
            // result = {u'A': 1, u'C': 4, u'$': 0, u'G': 7, u'T': 8}

            // PY: last_column, suffix_array = BWT(Text)
            // sorted(range(SZ), key=lambda p: text[p:]) -> indeksy posortowanych kawalkow
            // bwt = [(text[(i - 1 + SZ) % SZ], i) for i in sorted(range(SZ), key=lambda p: text[p:])]
            // bwt = (znak poprzedzajacy indeks i, indeks i)
            // bwt = [(u'T', 12), (u'$', 0), (u'T', 7), (u'C', 2), (u'A', 1), (u'G', 5), (u'A', 8), (u'T', 4), (u'T', 11), (u'C', 6), (u'A', 3), (u'T', 10), (u'C', 9)]

            // bwt, sufarr = map(list, zip(*bwt))
            // rozzipowuje na dwie sekwencje:
            // bwt = [u'T', u'$', u'T', u'C', u'A', u'G', u'A', u'T', u'T', u'C', u'A', u'T', u'C']
            // sufarr = [12, 0, 7, 2, 1, 5, 8, 4, 11, 6, 3, 10, 9]

            // PY: count = Count(last_column)
            // zwraca slownik, dla kazdego symbolu {$, A, C, G, T}
            // '$': [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            // 'A': [0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3] ...
            // number w tablicy jest inkrementowany po kazdym pojawieniu sie symbolu w poprzedniej pozycji.
            // T jest pierwsze w last_column
            // 'T' : [0, 1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 5, 5]
            // rozmiar sekwencji jest +1 ponad dlugosc last_column (pierwszy element z definicji = 0)

            // potrzebne sa:
            // first occurence => tablica 5 elementow, index w posortowanym cyclic, przez zliczenie
            // last_column => tablica dlugosci tekstu
            // === do skrocenia
            // count => slownik, 5 kluczy, 5 tablic o dlugosci last_column + 1
            // suffix_array => tablica dlugosci tekstu
#if 0
            0  1  2  3  4  5  6  7  8  9 10 11 12
text        A  C  A  T  G  C  T  A  C  T  T  T  $

first col   $  A  A  A  C  C  C  G  T  T  T  T  T

suffarr    12  0  7  2  1  5  8  4 11  6  3 10  9  (index w text)
last col    T  $  T  C  A  G  A  T  T  C  A  T  C

count($) 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
count(A) 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3
count(C) 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3
count(G) 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1
count(T) 0, 1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 5, 5

first occ   {u'A': 1, u'C': 4, u'$': 0, u'G': 7, u'T': 8}
#endif

            m_subchroma.push_back('$');
            const auto & subchroma = m_subchroma;

            const auto full_suffix_array = BW::full_suffix_array(std::string("ACATGCTACTTT$"));
            const auto last_col = BW::last_column(std::string("ACATGCTACTTT$"), full_suffix_array);
            const auto count = BW::count(last_col, 1);
            const auto count3 = BW::count(last_col, 3);

            std::vector<BW::count_type> foo(last_col.size() + 1);
            for (std::size_t ix{0}; ix < foo.size(); ++ix)
            {
                foo[ix] = count3.value(4, ix, last_col);
            }
            const auto first_occ = BW::first_occurences(count3, last_col);
//            const auto full_suffix_array = BW::full_suffix_array(subchroma);
            ;
            // clear it prior to collecting another one
            std::cout << m_subchroma.size() << std::endl;
        }

        m_subchroma.clear();
    }

    int passReferenceGenome(
        int chromatidSequenceId,
        const std::vector<std::string> & chromatidSequence)
    {
        m_curr_chromaid = chromatidSequenceId;

        for (const auto & chroma_piece : chromatidSequence)
        {
            std::vector<std::string> Npieces;
            boost::split(Npieces, chroma_piece, [](char c){return c == 'N';}, boost::token_compress_on);
            if (Npieces.back().empty())
            {
                Npieces.pop_back();
            }
            if (Npieces.front().empty())
            {
                Npieces.erase(Npieces.begin());
            }


            if (chroma_piece.front() == 'N')
            {
                close_sub_chroma();
            }

            if (Npieces.empty() == false)
            {
                // process all excluding the last one
                for (std::size_t ix{0}; ix < Npieces.size() - 1; ++ix)
                {
                    m_subchroma.append(Npieces[ix]);
                    close_sub_chroma();
                }

                // do the last one, but without closing
                m_subchroma.append(Npieces.back());
            }

            if (chroma_piece.back() == 'N')
            {
                close_sub_chroma();
            }
        }

        return 0;
    }

    std::vector<std::string> getAlignment(
        int N,
        double normA,
        double normS,
        const std::vector<std::string> & readName,
        const std::vector<std::string> & readSequence)
    {
        std::vector<std::string> ret(N, "");
        for (int i = 0; i < N; ++i)
        {
            std::string qname = "sim" + std::to_string(1 + i / 2) + '/'
                + ((i % 2) ? '2' : '1');
            ret[i] = qname + ",20,1,150,+,0";
        }
        return ret;
    }

    int m_curr_chromaid;
    std::string m_subchroma;
};

#endif /* DNA_SEQUENCING_HPP_ */
