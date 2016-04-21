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
#include <valarray>


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


double norm_pdf(double x, double mu, double sigma)
{
    const double arg = (x - mu) / sigma;
    return std::exp(- arg * arg / 2) / (sigma * std::sqrt(2. * M_PI));
}

double score_close_pairs(const std::vector<std::tuple<int, int, int>> & close_pairs)
{
    std::valarray<double> probs(close_pairs.size());

    std::transform(close_pairs.cbegin(), close_pairs.cend(), std::begin(probs),
        [](const std::tuple<int, int, int> & rpair)
        {
            const auto prob = norm_pdf(std::abs(std::get<1>(rpair) - std::get<2>(rpair)), 450., 34.);
            return prob;
        });

    return probs[0] / probs.sum();
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
            if (ix % 2000 == 0)
            {
                std::cerr << "[DNAS1] Doing read pair " << ix / 2 + 1 << " out of " << readName.size() / 2 << std::endl;
            }

            const auto & head_name = readName[ix];
            const auto & tail_name = readName[ix + 1];
            const auto & head_read_fwd = readSequence[ix];
            const auto & tail_read_fwd = readSequence[ix + 1];
            const auto head_read_rev = reverse_complement(head_read_fwd);
            const auto tail_read_rev = reverse_complement(tail_read_fwd);

            //              <chroma_id, head_pos, tail_pos>
            // range of distances between true positives: [263,814]
            // distance = (tail_center - head_center)
            enum
            {
                REAL_MIN_DIST = 263,
                REAL_MAX_DIST = 814,
                MIN_DIST = 450 - 300,
                MAX_DIST = 450 + 300,
            };
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
                        if (dist >= MIN_DIST && dist <= REAL_MAX_DIST)
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
                        if (dist >= MIN_DIST && dist <= REAL_MAX_DIST)
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

                const std::string confidence = std::to_string(score_close_pairs(close_pairs));
//                const std::string confidence = std::to_string(1. / close_pairs.size());

                head_res = head_name + ',' + std::to_string(std::get<0>(close_pairs.front())) + ',' +
                    std::to_string(std::abs(std::get<1>(close_pairs.front())) + 1) + ',' +
                    std::to_string(std::abs(std::get<1>(close_pairs.front())) + 150) + ',' +
                    (std::get<1>(close_pairs.front()) > 0 ? '+' : '-') + ',' + confidence;

                tail_res = tail_name + ',' + std::to_string(std::get<0>(close_pairs.front())) + ',' +
                    std::to_string(std::abs(std::get<2>(close_pairs.front())) + 1) + ',' +
                    std::to_string(std::abs(std::get<2>(close_pairs.front())) + 150) + ',' +
                    (std::get<2>(close_pairs.front()) > 0 ? '-' : '+') + ',' + confidence;
            }
            else if (cumm_matched_head_fwd.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_head_fwd.front());
                const BW::pos_type position = std::get<1>(cumm_matched_head_fwd.front());

                const std::string confidence_h = "0.97"; // approx 0.9741749322548726
                const std::string confidence_t = "0.49";
//                const std::string confidence_h = std::to_string(1. / cumm_matched_head_fwd.size());
//                const std::string confidence_t = std::to_string(1. / cumm_matched_head_fwd.size());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '+' : '-') + ',' + confidence_h;

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 451) + ',' +
                    std::to_string(position + 600) + ',' +
                    (position > 0 ? '-' : '+') + ',' + confidence_t;
            }
            else if (cumm_matched_tail_rev.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_tail_rev.front());
                const BW::pos_type position = std::get<1>(cumm_matched_tail_rev.front());

                const std::string confidence_h = "0.49";
                const std::string confidence_t = "0.97";
//                const std::string confidence_h = std::to_string(1. / cumm_matched_head_fwd.size());
//                const std::string confidence_t = std::to_string(1. / cumm_matched_head_fwd.size());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 - 450) + ',' +
                    std::to_string(position + 150 - 450) + ',' +
                    (position > 0 ? '+' : '-') + ',' + confidence_h;

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '-' : '+') + ',' + confidence_t;
            }
            else if (cumm_matched_head_rev.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_head_rev.front());
                const BW::pos_type position = std::get<1>(cumm_matched_head_rev.front());

                const std::string confidence_h = "0.97";
                const std::string confidence_t = "0.49";
//                const std::string confidence_h = std::to_string(1. / cumm_matched_head_fwd.size());
//                const std::string confidence_t = std::to_string(1. / cumm_matched_head_fwd.size());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '-' : '+') + ',' + confidence_h;

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 - 450) + ',' +
                    std::to_string(position - 300) + ',' +
                    (position > 0 ? '+' : '-') + ',' + confidence_t;
            }
            else if (cumm_matched_tail_fwd.size())
            {
                const chrid_type chroma_id = std::get<0>(cumm_matched_tail_fwd.front());
                const BW::pos_type position = std::get<1>(cumm_matched_tail_fwd.front());

                const std::string confidence_h = "0.49";
                const std::string confidence_t = "0.97";
//                const std::string confidence_h = std::to_string(1. / cumm_matched_head_fwd.size());
//                const std::string confidence_t = std::to_string(1. / cumm_matched_head_fwd.size());

                head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1 + 450) + ',' +
                    std::to_string(position + 600) + ',' +
                    (position > 0 ? '-' : '+') + ',' + confidence_h;

                tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                    std::to_string(position + 1) + ',' +
                    std::to_string(position + 150) + ',' +
                    (position > 0 ? '+' : '-') + ',' + confidence_t;
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
