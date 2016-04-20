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

#include "BW/BurrowsWheeler.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <utility>
#include <unordered_map>
#include <cstddef>


#include <sys/time.h>
long int timestamp()
{
    timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec;// + tv.tv_usec / 1e6;
}


std::string reverse_complement(const std::string & seq)
{
    std::string result;
    result.reserve(seq.size());

    std::transform(seq.crbegin(), seq.crend(), std::back_inserter(result),
        [](const char base)
        {
            return BW::nucleobase_by_ix(BW::ix_by_nucleobase(base) ^ 0x3);
        }
    );

    return result;
}


struct DNASequencing
{
    enum { SKIP = 100 };

    typedef int chrid_type;

    void initTest()
    {
        m_bw_contexts.clear();
    }

    int initTest(int testDifficulty)
    {
        initTest();
        return 0;
    }

    int preProcessing()
    {
        return 0;
    }

//#pragma GCC optimize ( "-O0" )

    int passReferenceGenome(
        int chromatidSequenceId,
        const std::vector<std::string> & chromatidSequence)
    {
        const auto time0 = timestamp();

        std::string chromatid;
        chromatid.reserve(chromatidSequence.size() * 100 + 1);

        for (const auto & chroma_piece : chromatidSequence)
        {
            std::string acgt_chroma = chroma_piece;

            std::transform(
                chroma_piece.cbegin(), chroma_piece.cend(),
                acgt_chroma.begin(),
                [](char const what)
                {
                    if (what != 'A' && what != 'C' && what != 'G' && what != 'T')
                    {
                        return 'N';
                    }
                    else
                    {
                        return what;
                    }
                });
            chromatid.append(acgt_chroma);
        }
        std::cerr << "[DNAS1] Combined chromatid length " << chromatid.size() << std::endl;
        chromatid.push_back('$');

        auto full_suffix_array = BW::full_suffix_array(chromatid);
        std::cerr << "[DNAS1] suffix array. Done " << std::endl;
        const auto last_col = BW::last_column(chromatid, full_suffix_array);
        const auto count = BW::count(last_col, SKIP);
        const auto first_occ = BW::first_occurences(count, last_col);
        const auto partial_sufarr = BW::partial_suffix_array(full_suffix_array, SKIP);
        full_suffix_array.clear();


        m_bw_contexts.emplace(chromatidSequenceId, BW::Context{last_col, count, first_occ, partial_sufarr});

        std::cerr << "[DNAS1] elapsed time (passReferenceGenome) " << timestamp() - time0 << " secs" << std::endl;

        return 0;
    }


//#pragma GCC optimize ( "-O0" )

    std::vector<std::string> getAlignment(
        int N,
        double normA,
        double normS,
        const std::vector<std::string> & readName,
        const std::vector<std::string> & readSequence)
    {
        const auto time0 = timestamp();

        std::vector<std::string> ret;
        ret.reserve(N);

        for (std::size_t ix{0}; ix < readName.size(); ix += 2)
        {
            if (ix % 100 == 0)
            {
                std::cerr << "Doing read pair " << ix / 2 + 1 << " out of " << readName.size() / 2 << std::endl;
            }

            const auto & head_name = readName[ix];
            const auto & tail_name = readName[ix + 1];
            const auto & head_read_fwd = readSequence[ix];
            const auto & tail_read_fwd = readSequence[ix + 1];
            const auto head_read_rev = reverse_complement(head_read_fwd);
            const auto tail_read_rev = reverse_complement(tail_read_fwd);

            //              <chroma_id, head_pos, tail_pos>
            std::vector<std::tuple<chrid_type, int, int>> close_pairs;

            //             <chroma_id, position>
            std::vector<std::pair<chrid_type, BW::pos_type>> cumm_matched_head_fwd;
            std::vector<std::pair<chrid_type, BW::pos_type>> cumm_matched_head_rev;
            std::vector<std::pair<chrid_type, BW::pos_type>> cumm_matched_tail_fwd;
            std::vector<std::pair<chrid_type, BW::pos_type>> cumm_matched_tail_rev;

            for (const auto & id_bw_context : m_bw_contexts)
            {
                const auto chroma_id = id_bw_context.first;
                const auto & bw_ctx = id_bw_context.second;

                const auto matched_head_fwd = BW::better_match(head_read_fwd, bw_ctx);
                const auto matched_head_rev = BW::better_match(head_read_rev, bw_ctx);
                const auto matched_tail_fwd = BW::better_match(tail_read_fwd, bw_ctx);
                const auto matched_tail_rev = BW::better_match(tail_read_rev, bw_ctx);

                for (const auto h_pos : matched_head_fwd)
                {
                    for (const auto t_pos : matched_tail_rev)
                    {
                        const auto dist = t_pos - h_pos;
                        if (dist >= 150 && dist <= 750)
                        {
                            close_pairs.emplace_back(chroma_id, h_pos, t_pos);
                        }
                    }
                }

                for (const auto h_pos : matched_head_rev)
                {
                    for (const auto t_pos : matched_tail_fwd)
                    {
                        const auto dist = h_pos - t_pos;
                        if (dist >= 150 && dist <= 750)
                        {
                            close_pairs.emplace_back(chroma_id, -h_pos, -t_pos);
                        }
                    }
                }

                auto append = [](std::vector<std::pair<chrid_type, BW::pos_type>> & cumm,
                                 const std::vector<BW::pos_type> & src,
                                 chrid_type chroma_id)
                {
                    cumm.reserve(cumm.size() + src.size());
                    std::transform(src.cbegin(), src.cend(), std::back_inserter(cumm),
                        [&chroma_id](const BW::pos_type what) -> std::pair<chrid_type, BW::pos_type>
                        {
                            return {chroma_id, what};
                        }
                    );
                };
                append(cumm_matched_head_fwd, matched_head_fwd, chroma_id);
                append(cumm_matched_head_rev, matched_head_rev, chroma_id);
                append(cumm_matched_tail_fwd, matched_tail_fwd, chroma_id);
                append(cumm_matched_tail_rev, matched_tail_rev, chroma_id);
            }


            std::string head_res;
            std::string tail_res;

            if (close_pairs.size() != 0)
            {
                std::sort(close_pairs.begin(), close_pairs.end(),
                    [](const std::tuple<chrid_type, int, int> & lhs, const std::tuple<chrid_type, int, int> & rhs)
                    {
                        const auto dlhs = std::abs(std::get<1>(lhs) - std::get<2>(lhs));
                        const auto drhs = std::abs(std::get<1>(rhs) - std::get<2>(rhs));
                        return std::abs(dlhs - 450) < std::abs(drhs - 450);
                    });

                head_res = head_name + ',' + std::to_string(std::get<0>(close_pairs.front())) + ',' +
                    std::to_string(std::abs(std::get<1>(close_pairs.front())) + 1) + ',' +
                    std::to_string(std::abs(std::get<1>(close_pairs.front())) + 150) + ',' +
                    (std::get<1>(close_pairs.front()) > 0 ? '+' : '-') + ",1.00";

                tail_res = tail_name + ',' + std::to_string(std::get<0>(close_pairs.front())) + ',' +
                    std::to_string(std::abs(std::get<2>(close_pairs.front())) + 1) + ',' +
                    std::to_string(std::abs(std::get<2>(close_pairs.front())) + 150) + ',' +
                    (std::get<2>(close_pairs.front()) > 0 ? '-' : '+') + ",1.00";
            }
            else if (cumm_matched_head_fwd.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_head_fwd.front());
                const BW::pos_type position = std::get<1>(cumm_matched_head_fwd.front());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '+' : '-') + ",0.99";

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 451) + ',' +
                    std::to_string(position + 600) + ',' +
                    (position > 0 ? '-' : '+') + ",0.50";
            }
            else if (cumm_matched_tail_rev.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_tail_rev.front());
                const BW::pos_type position = std::get<1>(cumm_matched_tail_rev.front());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 - 450) + ',' +
                    std::to_string(position + 150 - 450) + ',' +
                    (position > 0 ? '+' : '-') + ",0.50";

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '-' : '+') + ",0.99";
            }
            else if (cumm_matched_head_rev.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_head_rev.front());
                const BW::pos_type position = std::get<1>(cumm_matched_head_rev.front());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '-' : '+') + ",0.99";

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 - 450) + ',' +
                    std::to_string(position - 300) + ',' +
                    (position > 0 ? '+' : '-') + ",0.50";
            }
            else if (cumm_matched_tail_fwd.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_tail_fwd.front());
                const BW::pos_type position = std::get<1>(cumm_matched_tail_fwd.front());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 + 450) + ',' +
                    std::to_string(position + 600) + ',' +
                    (position > 0 ? '-' : '+') + ",0.50";

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '+' : '-') + ",0.99";
            }
            else
            {
                head_res = head_name + ",20,1,150,+,0.00";
                tail_res = tail_name + ",20,450,600,-,0.00";
            }

            ret.push_back(head_res);
            ret.push_back(tail_res);
        }

        std::cerr << "[DNAS1] elapsed time (getAlignment) " << timestamp() - time0 << " secs" << std::endl;

        return ret;
    }

    std::unordered_map<chrid_type, BW::Context> m_bw_contexts;
};

#endif /* DNA_SEQUENCING_HPP_ */
