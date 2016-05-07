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

//#pragma GCC optimize ( "-O3" )

#ifndef DNA_SEQUENCING_HPP_
#define DNA_SEQUENCING_HPP_

#define STORE_COMPRESSED_TEXT


#include "BW/BurrowsWheeler.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <cstddef>
#include <valarray>
#include <cmath>
#include <memory>

#include <stdlib.h>

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


typedef int chrid_type;

std::size_t LevenshteinDistance(const std::string& s1, const std::string& s2, const std::size_t max_dist = 150)
{
    typedef std::size_t size_type;

    const char * str1 = s1.c_str();
    const size_type SZ1 = s1.size();
    const char * str2 = s2.c_str();
    const size_type SZ2 = s2.size();


    typedef std::uint8_t dist_type;

    dist_type dist1[std::max(SZ1, SZ2) + 1];
    dist_type dist2[std::max(SZ1, SZ2) + 1];

    dist_type * v0{dist1};
    dist_type * v1{dist2};

    std::iota(v0, v0 + sizeof (dist1) / sizeof (dist1[0]), 0);

    for (size_type P{0}; P < SZ1; ++P)
    {
        v1[0] = P + 1;
        v1[1] = P + 2;

        dist_type mid_dist = max_dist;

        for (size_type Q{0}; Q < SZ2; ++Q)
        {
            const dist_type DIAG_SCORE = v0[Q] + (str1[P] != str2[Q]);

            const dist_type UP_DOWN_SCORE = std::min(v1[Q] + 1, v0[Q + 1] + 1);

            v1[Q + 1] = std::min(UP_DOWN_SCORE, DIAG_SCORE);
            mid_dist = std::min(v1[Q + 1], mid_dist);
        }

        std::swap(v0, v1);

        if (mid_dist > max_dist)
        {
            return max_dist + 1;
        }
    }

    const size_type result = v0[SZ2];

    return result;
}




std::size_t chromatid_offset(int test_type, int chroma_id)
{
    assert(test_type >= 0 && test_type < 3);
    assert(chroma_id >= 1 && chroma_id <= 24);

    --chroma_id;

    static const std::size_t offsets[][24] =
    {
        {
            // 20
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
        },
        {
            // 1, 11, 20
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            248956422u, 0, 0, 0, 0, 0, 0, 0, 0, 248956422u + 135086622u,
            0, 0, 0, 0,
        },
        {
            0u,
            248956422u,
            491149951u,
            689445510u,
            879660065u,
            1061198324u,
            1232004303u,
            1391350276u,
            1536488912u,
            1674883629u,
            1808681051u,
            1943767673u,
            2077042982u,
            2191407310u,
            2298451028u,
            2400442217u,
            2490780562u,
            2574038003u,
            2654411288u,
            2713028904u,
            2777473071u,
            2824183054u,
            2875001522u,
            3031042417u,
        }
    };

    return offsets[test_type][chroma_id];
}


std::pair<int, std::size_t> chroma_id_by_offset(std::size_t offset, int test_type)
{
    static const std::size_t offsets[][1 + 24] =
    {
        {
            0,
            // 20
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 64444167u,
            64444167u, 64444167u, 64444167u, 64444167u,
        },
        {
            0,
            // 1, 11, 20
            248956422u, 248956422u, 248956422u, 248956422u, 248956422u,
            248956422u, 248956422u, 248956422u, 248956422u, 248956422u,

            248956422u + 135086622u, 248956422u + 135086622u,
            248956422u + 135086622u, 248956422u + 135086622u,
            248956422u + 135086622u,
            248956422u + 135086622u, 248956422u + 135086622u,
            248956422u + 135086622u, 248956422u + 135086622u,
            248956422u + 135086622u + 64444167u,
            248956422u + 135086622u + 64444167u, 248956422u + 135086622u + 64444167u,
            248956422u + 135086622u + 64444167u, 248956422u + 135086622u + 64444167u,
        },
        {
            0,
            248956422u,
            491149951u,
            689445510u,
            879660065u,
            1061198324u,
            1232004303u,
            1391350276u,
            1536488912u,
            1674883629u,
            1808681051u,
            1943767673u,
            2077042982u,
            2191407310u,
            2298451028u,
            2400442217u,
            2490780562u,
            2574038003u,
            2654411288u,
            2713028904u,
            2777473071u,
            2824183054u,
            2875001522u,
            3031042417u,
            3031042417u + 57227415u
        }
    };

    const auto & arr = offsets[test_type];

    for (std::size_t ix = 1; ix <= 24; ++ix)
    {
        if (offset <= arr[ix])
        {
            return {ix, arr[ix - 1]};
        }
    }

    return {24, 3031042417u};
}


std::size_t genome_reserved_space(int const test_type)
{
#if 0
1   248956422
2   242193529
3   198295559
4   190214555
5   181538259
6   170805979
7   159345973
8   145138636
9   138394717
10  133797422
11  135086622
12  133275309
13  114364328
14  107043718
15  101991189
16  90338345
17  83257441
18  80373285
19  58617616
20  64444167
21  46709983
22  50818468
23  156040895
24  57227415
#endif
    static const std::size_t reserved[] =
    {
        1ULL + 64444167,
        1ULL + 248956422 + 135086622 + 64444167,
        1ULL +
        248956422 +
        242193529 +
        198295559 +
        190214555 +
        181538259 +
        170805979 +
        159345973 +
        145138636 +
        138394717 +
        133797422 +
        135086622 +
        133275309 +
        114364328 +
        107043718 +
        101991189 +
        90338345 +
        83257441 +
        80373285 +
        58617616 +
        64444167 +
        46709983 +
        50818468 +
        156040895 +
        57227415
    };
    constexpr long RESERVED_SZ = sizeof (reserved) / sizeof (reserved[0]);

    const std::size_t retval =
        (test_type >= 0 && test_type < RESERVED_SZ) ? reserved[test_type] : 0;

    return retval;
}


struct DNASequencing
{
    enum { SKIP = 32 };

    int m_test_type;

    std::string m_genome;

    BW::Context m_bw_context;

#ifdef STORE_COMPRESSED_TEXT
    BW::CompressedText m_compressed_text;
#endif


    // 1.
    int initTest(int testDifficulty)
    {
        m_test_type = testDifficulty;

        m_genome.resize(genome_reserved_space(m_test_type));

        return 0;
    }


    // 2.
    int passReferenceGenome(
        int chromatidSequenceId,
        const std::vector<std::string> & chromatidSequence);


    // 3.
    int preProcessing();


    // 4.
    std::vector<std::string> getAlignment(
        int N,
        double normA,
        double normS,
        const std::vector<std::string> & readName,
        const std::vector<std::string> & readSequence);
};





int
DNASequencing::passReferenceGenome(
    int chromatidSequenceId,
    const std::vector<std::string> & chromatidSequence)
{
    std::cerr << "[DNAS1] passReferenceGenome start" << std::endl;
    const auto time0 = timestamp();

    std::size_t local_offset{0};

    for (const auto & kmer : chromatidSequence)
    {
        std::transform(
            kmer.cbegin(), kmer.cend(),
            m_genome.begin() + chromatid_offset(m_test_type, chromatidSequenceId) + local_offset,
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

        local_offset += kmer.size();
    }

    std::cerr << "[DNAS1] elapsed time (passReferenceGenome) " << timestamp() - time0 << " secs" << std::endl;

    return 0;
}


int DNASequencing::preProcessing()
{
    const auto time0 = timestamp();

    m_genome.back() = '$';

    std::cerr << "[DNAS1] Combined genome length " << m_genome.size() << std::endl;

#ifdef STORE_COMPRESSED_TEXT
    m_compressed_text = BW::compress_text(m_genome);
    std::cerr << "[DNAS1] compressed text. Done" << std::endl;
#endif

    auto full_suffix_array = BW::full_suffix_array(m_genome);
    std::cerr << "[DNAS1] suffix array. Done" << std::endl;
    const auto last_col = BW::last_column(m_genome, full_suffix_array);
    std::cerr << "[DNAS1] last column. Done" << std::endl;
    m_genome.clear();
    const auto partial_sufarr = BW::partial_suffix_array(full_suffix_array, SKIP);
    std::cerr << "[DNAS1] partial suffix array. Done" << std::endl;
    full_suffix_array.clear();

    const auto count = BW::count(last_col, SKIP);
    const auto first_occ = BW::first_occurences(count, last_col);
    m_bw_context = BW::Context{last_col, count, first_occ, partial_sufarr};

    std::cerr << "[DNAS1] elapsed time (preProcessing) " << timestamp() - time0 << " secs" << std::endl;

    return 0;
}


int64_t memcmpr(const char * a, const char *b, std::size_t n)
__attribute__ ((__pure__))
__attribute__ ((__nonnull__ (1, 2)));

int64_t memcmpr(const char * a, const char * b, std::size_t n)
{
    auto tail = n % 8;
    const auto head = n - tail;
    while (tail--)
    {
        if (a[head + tail] != b[head + tail])
        {
            return a[head + tail] - b[head + tail];
        }
    }

    const int64_t * a8 = (const int64_t *)a;
    const int64_t * b8 = (const int64_t *)b;
    n = n / 8;

    while (n--)
    {
        if (a8[n] != b8[n])
        {
            return a8[n] - b8[n];
        }
    }

    return 0;
}


std::vector<std::string>
DNASequencing::getAlignment(
    int N,
    double normA,
    double normS,
    const std::vector<std::string> & readName,
    const std::vector<std::string> & readSequence)
{
    const auto time0 = timestamp();

    // create reverse-complements of reads
    enum {READ_SZ = 150};
    const std::vector<std::string> & reads_fwd{readSequence};
    std::vector<std::string> reads_rev;

    reads_rev.reserve(reads_fwd.size());
    std::transform(reads_fwd.cbegin(), reads_fwd.cend(), std::back_inserter(reads_rev),
        reverse_complement);


    ////////////////////////////////////////////////////////////////////////////


    // create kmer views
    enum {KMER_SZ = 75};

    static_assert((READ_SZ % KMER_SZ) == 0, "KMER_SZ doesn't divide READ_SZ");
    enum {KMER_PER_READ = READ_SZ / KMER_SZ};

    typedef uint32_t rkix_type;
    const rkix_type NEG_RKIX_OFFSET = readName.size() * KMER_PER_READ;

    const auto make_rkix = [NEG_RKIX_OFFSET](uint32_t rix, uint8_t kmer_pos, bool is_reverse) -> rkix_type
    {
        return (rix * READ_SZ + kmer_pos) / KMER_SZ + (is_reverse ? NEG_RKIX_OFFSET : 0);
    };

    std::vector<rkix_type> kmer_views;

    // twice (fwd + rev) number of kmers per read
    kmer_views.reserve(2 * KMER_PER_READ * readName.size());

    for (std::size_t rix = 0; rix < reads_fwd.size(); ++rix)
    {
        for (std::size_t kix = 0; kix < READ_SZ; kix += KMER_SZ)
        {
            const auto rkix = make_rkix(rix, kix, false);
            kmer_views.push_back(rkix);
        }
    }
    for (std::size_t rix = 0; rix < reads_rev.size(); ++rix)
    {
        for (std::size_t kix = 0; kix < READ_SZ; kix += KMER_SZ)
        {
            const auto rkix = make_rkix(rix, kix, true);
            kmer_views.push_back(rkix);
        }
    }
    std::cerr << "[DNAS1] processing " << kmer_views.size() << " K-mer views" << std::endl;
    std::cerr << "[DNAS1] reverse-complement kmers, elapsed " << timestamp() - time0 << " secs" << std::endl;


    ////////////////////////////////////////////////////////////////////////////


    // sort-em
    struct sort_ctx_type
    {
        const std::string * fwd;
        const std::string * rev;
        const rkix_type neg_rkix_offset;
    } sort_ctx {reads_fwd.data(), reads_rev.data(), NEG_RKIX_OFFSET};

    qsort_r(kmer_views.data(), kmer_views.size(), sizeof (kmer_views.front()),
        [](const void *va, const void * vb, void * vc)
        {
            const rkix_type * a = (const rkix_type *)va;
            const rkix_type * b = (const rkix_type *)vb;
            const sort_ctx_type * ctx = (const sort_ctx_type *)vc;

            const auto rkix_a = *a >= ctx->neg_rkix_offset ? *a - ctx->neg_rkix_offset : *a;
            const auto rix_a = (rkix_a * KMER_SZ) / READ_SZ;
            const auto kix_a = (rkix_a * KMER_SZ) % READ_SZ;

            const char * lhs = *a >= ctx->neg_rkix_offset ?
                ctx->rev[rix_a].c_str() + kix_a :
                ctx->fwd[rix_a].c_str() + kix_a;

            const auto rkix_b = *b >= ctx->neg_rkix_offset ? *b - ctx->neg_rkix_offset : *b;
            const auto rix_b = (rkix_b * KMER_SZ) / READ_SZ;
            const auto kix_b = (rkix_b * KMER_SZ) % READ_SZ;

            const char * rhs = *b >= ctx->neg_rkix_offset ?
                ctx->rev[rix_b].c_str() + kix_b :
                ctx->fwd[rix_b].c_str() + kix_b;

            const auto result = memcmpr(lhs, rhs, KMER_SZ);

            return result > 0 ? 1 : result < 0 ? -1 : 0;
        }, &sort_ctx);

    std::cerr << "[DNAS1] sorted kmer views, elapsed " << timestamp() - time0 << " secs" << std::endl;


    ////////////////////////////////////////////////////////////////////////////


    std::vector<std::pair<BW::pos_type, BW::pos_type>> top_bottoms(kmer_views.size(), {0, m_bw_context.last_column.size() - 1});


    for (std::size_t depth = 0; depth < KMER_SZ; ++depth)
    {
        const auto tbix = kmer_views[0];
        const auto rkix = tbix >= NEG_RKIX_OFFSET ? tbix - NEG_RKIX_OFFSET : tbix;
        const auto rix = (rkix * KMER_SZ) / READ_SZ;
        const auto kix = (rkix * KMER_SZ) % READ_SZ;

        const char * kmer_p = tbix >= NEG_RKIX_OFFSET ?
                reads_rev[rix].c_str() + kix :
                reads_fwd[rix].c_str() + kix;

        auto prev_symbol = kmer_p[KMER_SZ - depth - 1];
        auto prev_top_bottom = top_bottoms[tbix];

        top_bottoms[tbix] = BW::top_botom_iter(prev_symbol, m_bw_context, prev_top_bottom.first, prev_top_bottom.second);

        for (std::size_t vix = 1; vix < kmer_views.size(); ++vix)
        {
            const auto tbix = kmer_views[vix];
            const auto rkix = tbix >= NEG_RKIX_OFFSET ? tbix - NEG_RKIX_OFFSET : tbix;
            const auto rix = (rkix * KMER_SZ) / READ_SZ;
            const auto kix = (rkix * KMER_SZ) % READ_SZ;

            if (top_bottoms[tbix] == BW::TOP_BOTTOM_INVALID)
            {
                continue;
            }

            const char * kmer_p = tbix >= NEG_RKIX_OFFSET ?
                    reads_rev[rix].c_str() + kix :
                    reads_fwd[rix].c_str() + kix;

            const auto symbol = kmer_p[KMER_SZ - depth - 1];

            if (prev_top_bottom == top_bottoms[tbix] && prev_symbol == symbol)
            {
                top_bottoms[tbix] = top_bottoms[kmer_views[vix - 1]];
            }
            else
            {
                prev_top_bottom = top_bottoms[tbix];
                prev_symbol = symbol;

                top_bottoms[tbix] = BW::top_botom_iter(symbol, m_bw_context, top_bottoms[tbix].first, top_bottoms[tbix].second);
            }
        }
    }

    std::cerr << "[DNAS1] BW on kmers, elapsed " << timestamp() - time0 << " secs" << std::endl;


    ////////////////////////////////////////////////////////////////////////////


    std::vector<std::string> ret;
    ret.reserve(N);

    for (std::size_t ix{0}; ix < readName.size(); ix += 2)
    {
        if (ix % 10000 == 0)
        {
            std::cerr << "[DNAS1] Doing read pair " << ix / 2 + 1 << " out of " << readName.size() / 2 << std::endl;
        }

        const auto & head_name = readName[ix];
        const auto & tail_name = readName[ix + 1];
        const auto & head_read_fwd = reads_fwd[ix];
        const auto & tail_read_fwd = reads_fwd[ix + 1];
        const auto & head_read_rev = reads_rev[ix];
        const auto & tail_read_rev = reads_rev[ix + 1];

        enum
        {
            REAL_MIN_DIST = 263,
            REAL_MAX_DIST = 814,
            MIN_DIST = 450 - 300,
            MAX_DIST = 450 + 300,
        };
        std::vector<std::tuple<int32_t, int32_t, std::size_t>> close_pairs; // pos, pos, distance

        ////////////////////////////////////////////////////////////////////////


        const auto cached_approximate_bw_match = [&top_bottoms, this, &make_rkix](
            const std::size_t rix, bool is_reverse)
        {
            const auto & bw_context = this->m_bw_context;

            std::unordered_set<BW::pos_type> positions;

            for (std::size_t kix = 0; kix < READ_SZ; kix += KMER_SZ)
            {
                const auto rkix = make_rkix(rix, kix, is_reverse);
                const auto & top_bottom = top_bottoms.at(rkix);

                if (top_bottom == BW::TOP_BOTTOM_INVALID)
                {
                    continue;
                }

                for (auto lcix = top_bottom.first; lcix <= top_bottom.second; ++lcix)
                {
                    const auto position = bw_context.partial_suffix_array.value(
                        lcix,
                        bw_context.last_column,
                        bw_context.count,
                        bw_context.first_occurences);

                    positions.emplace(position - kix);
                }
            }

            return positions;
            //return std::vector<BW::pos_type>(positions.cbegin(), positions.cend());
        };
        ////////////////////////////////////////////////////////////////////////

        const auto matched_head_fwd = cached_approximate_bw_match(ix, false);
        const auto matched_tail_fwd = cached_approximate_bw_match(ix + 1, false);
        const auto matched_head_rev = cached_approximate_bw_match(ix, true);
        const auto matched_tail_rev = cached_approximate_bw_match(ix + 1, true);

        for (const auto h_pos : matched_head_fwd)
        {
            for (const auto t_pos : matched_tail_rev)
            {
                const auto dist = t_pos - h_pos;
                if (dist >= MIN_DIST && dist <= REAL_MAX_DIST)
                {

                    close_pairs.emplace_back(h_pos, t_pos, 0);
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
                    close_pairs.emplace_back(-h_pos, -t_pos, 0);
                }
            }
        }


        if (close_pairs.size() > 1 && close_pairs.size() <= 3)
        {
            for (auto & pair : close_pairs)
            {
                const auto h_pos = std::get<0>(pair);
                const auto t_pos = std::get<1>(pair);
                if (h_pos > 0)
                {
                    const auto h_text = m_compressed_text.decompress(h_pos, h_pos + READ_SZ);
                    const auto h_dist = LevenshteinDistance(h_text, head_read_fwd);

                    const auto t_text = m_compressed_text.decompress(t_pos, t_pos + READ_SZ);
                    const auto t_dist = LevenshteinDistance(t_text, tail_read_rev);

                    std::get<2>(pair) = h_dist + t_dist;
                }
                else
                {
                    const auto h_text = m_compressed_text.decompress(-h_pos, -h_pos + READ_SZ);
                    const auto h_dist = LevenshteinDistance(h_text, head_read_rev);

                    const auto t_text = m_compressed_text.decompress(-t_pos, -t_pos + READ_SZ);
                    const auto t_dist = LevenshteinDistance(t_text, tail_read_fwd);

                    std::get<2>(pair) = h_dist + t_dist;
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////


        std::string head_res;
        std::string tail_res;

        if (close_pairs.size() != 0)
        {
            std::sort(close_pairs.begin(), close_pairs.end(),
                [](const std::tuple<int32_t, int32_t, std::size_t> & lhs,
                   const std::tuple<int32_t, int32_t, std::size_t> & rhs)
                {
                    const auto ldelta = std::abs(std::get<0>(lhs) - std::get<1>(lhs));
                    const auto rdelta = std::abs(std::get<0>(rhs) - std::get<1>(rhs));
                    const auto ldist = std::get<2>(lhs);
                    const auto rdist = std::get<2>(rhs);

                    if (ldist != rdist)
                    {
                        return ldist < rdist;
                    }
                    else
                    {
                        return std::abs(ldelta - 450) < std::abs(rdelta - 450);
                    }
                });

            const auto & best = close_pairs.front();
            const auto score = std::to_string(1. / close_pairs.size());

            const auto chroma_id_shift = chroma_id_by_offset(std::abs(std::get<0>(best)), m_test_type);
            const auto chroma_id = chroma_id_shift.first;
            const auto chroma_shift = chroma_id_shift.second;

            head_res = head_name + ',' + std::to_string(chroma_id) + ',' +
                std::to_string(std::abs(std::get<0>(best)) - chroma_shift + 1) + ',' +
                std::to_string(std::abs(std::get<0>(best)) - chroma_shift + 150) + ',' +
                (std::get<0>(best) > 0 ? '+' : '-') + ',' + score;

            tail_res = tail_name + ',' + std::to_string(chroma_id) + ',' +
                std::to_string(std::abs(std::get<1>(best)) - chroma_shift + 1) + ',' +
                std::to_string(std::abs(std::get<1>(best)) - chroma_shift + 150) + ',' +
                (std::get<1>(best) > 0 ? '-' : '+') + ',' + score;
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


#endif /* DNA_SEQUENCING_HPP_ */
