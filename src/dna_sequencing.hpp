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
#include "boost/tokenizer.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>


struct DNASequencing
{
    void initTest()
    {
        m_curr_chromaid = -1;
        m_chroma_off = 0;
        m_sub_chroma_off = 0;
        m_subchroma.clear();
        m_chromatid_bw_contexts.clear();
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

        std::cerr << "[DNAS1] Received " << m_chromatid_bw_contexts.size() << " chromatid(s)" << std::endl;
        for (const auto & pair : m_chromatid_bw_contexts)
        {
            std::cerr << "[DNAS1] Chromatid " << pair.first << ": " << pair.second.size() << " subchromatid(s)" << std::endl;
        }

        return 0;
    }

    void close_sub_chroma()
    {
        constexpr std::size_t MIN_SUBCHROMA_SZ = 300;

        // process subchroma collected so far
        if (m_subchroma.size() >= MIN_SUBCHROMA_SZ)
        {
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

                  *        *        *        *
            0  1  2  3  4  5  6  7  8  9 10 11 12
           ---------------------------------------
suffarr    12  0  7  2  1  5  8  4 11  6  3 10  9  (index w text)
           ---------------------------------------
suffarr           7        5     ? 11       10     (index w text)
'G'                              ^

#endif

            m_subchroma.push_back('$');
            const auto & subchroma = m_subchroma;

//            {
//                const auto full_suffix_array = BW::full_suffix_array(std::string("ACATGCTACTTT$"));
//                const auto last_col = BW::last_column(std::string("ACATGCTACTTT$"), full_suffix_array);
//                const auto count = BW::count(last_col, 1);
//                const auto count3 = BW::count(last_col, 3);
//
//                std::vector<BW::count_type> foo(last_col.size() + 1);
//                for (std::size_t ix{0}; ix < foo.size(); ++ix)
//                {
//                    foo[ix] = count3.value(4, ix, last_col);
//                }
//                const auto first_occ = BW::first_occurences(count3, last_col);
//                const auto a = BW::ix_by_base('$');
//                const auto b = BW::ix_by_base('A');
//                const auto c = BW::ix_by_base('C');
//                const auto d = BW::ix_by_base('G');
//                const auto e = BW::ix_by_base('T');
//                const auto partial_sufarr = BW::partial_suffix_array(full_suffix_array, 3);
//                std::vector<std::size_t> bar(full_suffix_array.size());
//                for (std::size_t ix{0}; ix < bar.size(); ++ix)
//                {
//                    bar[ix] = partial_sufarr.value(ix, last_col, count3, first_occ);
//                }
//                const auto res = BW::better_match("CATG", BW::Context{last_col, count3, first_occ, partial_sufarr});
//            }

            auto full_suffix_array = BW::full_suffix_array(subchroma);
            const auto last_col = BW::last_column(subchroma, full_suffix_array);
            const auto count = BW::count(last_col, 256);
            const auto first_occ = BW::first_occurences(count, last_col);
            const auto partial_sufarr = BW::partial_suffix_array(full_suffix_array, 256);
            full_suffix_array.clear();

            m_chromatid_bw_contexts[m_curr_chromaid].emplace_back(
                m_sub_chroma_off, BW::Context{last_col, count, first_occ, partial_sufarr});

//            std::cerr << "[DNAS1] " << m_subchroma.size() << std::endl;
        }
        if (m_subchroma.size()) // temporary for reference logging
        {
            std::cerr << "[DNAS1] " << m_subchroma.size() << ", off: " << m_sub_chroma_off << std::endl;
        }

        m_subchroma.clear();
        m_sub_chroma_off = m_chroma_off;
    }

    int passReferenceGenome(
        int chromatidSequenceId,
        const std::vector<std::string> & chromatidSequence)
    {
        if (m_curr_chromaid != chromatidSequenceId)
        {
            // reset offset with the chromatid when we start receiving a new one
            m_chroma_off = 0;
            m_sub_chroma_off = 0;
        }
        m_curr_chromaid = chromatidSequenceId;

        for (const auto & chroma_piece : chromatidSequence)
        {
            // tokenizer
            typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
            boost::char_separator<char> sep("N", "", boost::drop_empty_tokens);
            tokenizer tokens(chroma_piece, sep);

            std::vector<std::pair<uint32_t, std::string>> Npieces;

            for (tokenizer::iterator tok_iter = tokens.begin();
                     tok_iter != tokens.end();
                     ++tok_iter)
            {
                const uint32_t offset = tok_iter.base() - chroma_piece.begin() - tok_iter->size();

                Npieces.emplace_back(offset, *tok_iter);
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
                    m_sub_chroma_off += Npieces[ix].first;
                    m_subchroma.append(Npieces[ix].second);
                    close_sub_chroma();
                }

                // do the last one, but without closing
                m_sub_chroma_off += Npieces.back().first;
                m_subchroma.append(Npieces.back().second);
            }

            m_chroma_off += chroma_piece.size();

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
        for (std::size_t ix{0}; ix < readName.size(); ix += 2)
        {
            const auto & head_name = readName[ix];
            const auto & tail_name = readName[ix + 1];
            const auto & head_read = readSequence[ix];
            const auto & tail_read = readSequence[ix + 1];

            for (const auto & chroma_bw : m_chromatid_bw_contexts)
            {
                const auto chroma_id = chroma_bw.first;

                for (const auto & offset_bw_ctx : chroma_bw.second)
                {
                    const auto offset = offset_bw_ctx.first;
                    const auto & bw_ctx = offset_bw_ctx.second;
                    const auto matched = BW::better_match(head_read, bw_ctx);
                    if (matched.size())
                    {
                        std::cerr << "[DNAS1] matched " << matched.size() << std::endl;
                        std::cerr << "[DNAS1] " << matched.front() + offset << std::endl;
                        std::cerr << "[DNAS1] " << head_read << std::endl;
//                        std::strncmp(head_read.c_str(), );
                    }
                }
            }
        }


        std::vector<std::string> ret(N, "");

        for (int i = 0; i < N; ++i)
        {
            std::string qname = "sim" + std::to_string(1 + i / 2) + '/'
                + ((i % 2) ? '2' : '1');
            ret[i] = qname + ",20,1,150,+,0";
        }

        return ret;
    }

    // awful state, ugh:

    // chromatid id, passed with each passReferenceGenome
    int m_curr_chromaid;
    // global position within chromatid being collected (across calls to passReferenceGenome)
    // irrelevant of splits
    // shall be reset for each new chromatid
    std::size_t m_chroma_off;
    // global position within chromatid being collected (across calls to passReferenceGenome)
    // irrelevant of splits
    // which points to the start of subchromatid being collected
    std::size_t m_sub_chroma_off;
    // subsection of chromatid being collected,
    // subsections are separated by Ns in a chromatid
    std::string m_subchroma;
    // lists of BW contexts for each chromatid collected
    std::unordered_map<int, std::vector<std::pair<uint32_t, BW::Context>>> m_chromatid_bw_contexts;
};

#endif /* DNA_SEQUENCING_HPP_ */
