/*******************************************************************************
 * Copyright (c) 2016 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: BurrowsWheeler.hpp
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


#ifndef BURROWSWHEELER_HPP_
#define BURROWSWHEELER_HPP_

#include <cstdint>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <iterator>
#include <cstring>
#include <memory>
#include <utility>
#include <unordered_map>
#include <cstdint>
#include <cassert>
#include <set>


namespace BW
{


typedef uint32_t pos_type;


template <typename SeqT>
std::vector<pos_type>
full_suffix_array(const SeqT & text)
{
    std::vector<pos_type> sarray(text.size());

    std::iota(sarray.begin(), sarray.end(), 0);

    struct sort_ctx_type
    {
        const char * base;
        std::size_t sz;
    } sort_ctx {text.data(), text.size()};

    qsort_r(sarray.data(), sarray.size(), sizeof (sarray[0]),
        [](const void *va, const void * vb, void * vc)
        {
            const pos_type * a = (const pos_type *)va;
            const pos_type * b = (const pos_type *)vb;
            const sort_ctx_type * ctx = (const sort_ctx_type *)vc;

            typename SeqT::const_pointer lhs = ctx->base + *a;
            typename SeqT::const_pointer rhs = ctx->base + *b;

            //const std::size_t how_much = ctx->sz - std::max(*a, *b);
            const std::size_t how_much = std::min(ctx->sz - std::max(*a, *b), 1000UL);

            const auto result = std::memcmp(lhs, rhs, how_much);

            return result > 0 ? 1 : result < 0 ? -1 : 0;
        }, &sort_ctx);

    return sarray;
}


template <typename SeqT>
SeqT last_column(const SeqT & text, const std::vector<pos_type> & full_suffix_array)
{
//    SeqT last_col(text.size(), {});
    SeqT last_col;
    last_col.resize(text.size());

    assert(full_suffix_array.size() == last_col.size());

    std::transform(full_suffix_array.cbegin(), full_suffix_array.cend(),
        last_col.begin(),
        [&text](const pos_type ix)
        {
            return ix == 0 ? '$' : text[ix - 1];
        });

    return last_col;
}


constexpr inline
char base_by_ix(const std::size_t ix)
{
    return ((
        ((uint64_t)'N' << 40) |
        ((uint64_t)'$' << 32) |
        ((uint64_t)'T' << 24) |
        ((uint64_t)'G' << 16) |
        ((uint64_t)'C' << 8) |
        ((uint64_t)'A' << 0) |
        0
        ) >> (ix * 8)) & 0xFF;
}


constexpr inline
char nucleobase_by_ix(const std::size_t ix)
{
    return ((
        ((uint32_t)'T' << 24) |
        ((uint32_t)'G' << 16) |
        ((uint32_t)'C' << 8) |
        ((uint32_t)'A' << 0) |
        0
        ) >> (ix * 8)) & 0xFF;
}


constexpr inline
std::size_t ix_by_nucleobase(const char base)
{
    return
        ((((0ULL << ('A' - 'A')) |
            (1ULL << ('C' - 'A')) |
            (2ULL << ('G' - 'A')) |
            (3ULL << ('T' - 'A'))) >> (base - 'A')) & 0x3);
}


constexpr inline
std::size_t ix_by_base(const char base)
{
    return
        base == '$' ?
            4
            : base == 'N' ?
                5
                :
                ix_by_nucleobase(base);
}


typedef uint32_t count_type;

struct Count
{
    //const
    std::size_t m_skip;
    //const
    std::vector<std::vector<count_type>> m_counts;

    template <typename SeqT>
    inline
    count_type value(
        const std::size_t base_ix,
        const std::size_t position,
        const SeqT & last_column) const
    {
        assert(position <= last_column.size());

        if (position == 0)
        {
            return 0;
        }

        const auto remaining = position % m_skip;
        const auto mod_pos = position - remaining;
        const auto mod_count = m_counts[base_ix][mod_pos / m_skip];

        if (remaining == 0)
        {
            return mod_count;
        }

        const auto base = base_by_ix(base_ix);

        const auto real_count =
            std::accumulate(
                std::next(last_column.cbegin(), mod_pos),
                std::next(last_column.cbegin(), position),
                mod_count,
                [&base](count_type acc, typename SeqT::value_type item)
                {
                    return acc + (item == base);
                }
            );

        return real_count;
    }
};


template <typename SeqT>
Count count(const SeqT & last_column, const std::size_t SKIP = 4)
{
    std::vector<std::vector<count_type>> result(6);
    std::vector<count_type> stat(6);

    for (auto & vec : result)
    {
        vec.push_back(0);
    }

    for (std::size_t ix{0}; ix < last_column.size(); ++ix)
    {
        const auto item = last_column[ix];

        for (std::size_t six{0}; six < stat.size(); ++six)
        {
            stat[six] += (item == base_by_ix(six));
        }

        if (ix % SKIP == (SKIP - 1))
        {
            for (std::size_t rix{0}; rix < result.size(); ++rix)
            {
                result[rix].push_back(stat[rix]);
            }
        }
    }

    return Count{SKIP, result};
}


template <typename SeqT>
std::vector<std::size_t>
first_occurences(const Count & count, const SeqT & last_column)
{
    std::vector<std::size_t> first_occ(6);

    first_occ[ix_by_base('$')] = 0;
    first_occ[ix_by_base('A')] = 1;
    first_occ[ix_by_base('C')] = count.value(ix_by_base('A'), last_column.size(), last_column) + first_occ[ix_by_base('A')];
    first_occ[ix_by_base('G')] = count.value(ix_by_base('C'), last_column.size(), last_column) + first_occ[ix_by_base('C')];
    first_occ[ix_by_base('N')] = count.value(ix_by_base('G'), last_column.size(), last_column) + first_occ[ix_by_base('G')];
    first_occ[ix_by_base('T')] = count.value(ix_by_base('N'), last_column.size(), last_column) + first_occ[ix_by_base('N')];

    return first_occ;
}


struct PartialSuffixArray
{
    //const
    std::size_t m_skip;
    //const
    std::unordered_map<pos_type, pos_type> m_sufarr;

    template <typename SeqT>
    pos_type value(
        pos_type ix,
        const SeqT & last_column,
        const Count & count,
        const std::vector<std::size_t> & first_occurrences) const
    {
        assert(ix < last_column.size());

        std::size_t backtrack{0};

        while (m_sufarr.count(ix) == 0)
        {
            assert(ix < last_column.size());
            const auto pred_base = last_column[ix];
            const auto base_ix = ix_by_base(pred_base);
            const auto base_ord = count.value(base_ix, ix, last_column);
            ix = first_occurrences[base_ix] + base_ord;
            ++backtrack;
        }

        return (m_sufarr.at(ix) + backtrack) % last_column.size();
    }
};

PartialSuffixArray
partial_suffix_array(const std::vector<pos_type> & full_suffix_array, const std::size_t SKIP = 4)
{
    if (full_suffix_array.size() == 0)
    {
        return PartialSuffixArray{SKIP, {}};
    }

    std::unordered_map<pos_type, pos_type> partial;
    for (std::size_t ix{0}; ix < full_suffix_array.size(); ix += SKIP)
    {
        partial[ix] = full_suffix_array[ix];
    }

    return PartialSuffixArray{SKIP, partial};
}


struct Context
{
    //const
    std::string last_column;
    //const
    BW::Count count;
    //const
    std::vector<std::size_t> first_occurences;
    //const
    BW::PartialSuffixArray partial_suffix_array;
};


static constexpr std::pair<pos_type, pos_type> TOP_BOTTOM_INVALID = {1, 0};


inline
std::pair<pos_type, pos_type>
top_botom_iter(
    const char symbol,
    const Context & ctx,
    const std::size_t top,
    const std::size_t bottom)
{
    const auto & last_column(ctx.last_column);
    const auto & count(ctx.count);
    const auto & first_occurences(ctx.first_occurences);

    // TODO memchr
    if (std::find(last_column.begin() + top, last_column.begin() + bottom + 1, symbol) != last_column.begin() + bottom + 1)
    {
        const auto base_ix = ix_by_base(symbol);

        const auto newtop = first_occurences.at(base_ix) + count.value(base_ix, top, last_column);
        const auto newbottom = first_occurences.at(base_ix) + count.value(base_ix, bottom + 1, last_column) - 1;

        return {newtop, newbottom};
    }
    else
    {
        return TOP_BOTTOM_INVALID;
    }
}


std::pair<pos_type, pos_type>
top_bottom(
    const char * begin,
    const char * end,
    const Context & ctx,
    std::size_t top,
    std::size_t bottom)
{
    const auto & last_column(ctx.last_column);
    const auto & count(ctx.count);
    const auto & first_occurences(ctx.first_occurences);

    auto pattern_ix{std::distance(begin, end)};
    assert(pattern_ix >= 0);

    while (top <= bottom)
    {
        if (pattern_ix != 0)
        {
            const auto symbol = begin[-1 + pattern_ix--];
            if (std::find(last_column.begin() + top, last_column.begin() + bottom + 1, symbol) != last_column.begin() + bottom + 1)
            {
                const auto base_ix = ix_by_base(symbol);

                top = first_occurences.at(base_ix) + count.value(base_ix, top, last_column);
                bottom = first_occurences.at(base_ix) + count.value(base_ix, bottom + 1, last_column) - 1;
            }
            else
            {
                return TOP_BOTTOM_INVALID;
            }
        }
        else
        {
            break;
        }
    }

    return {top, bottom};
}


std::pair<pos_type, pos_type>
top_bottom(
    const char * begin,
    const char * end,
    const Context & ctx)
{
    const auto & last_column(ctx.last_column);
    const auto & count(ctx.count);
    const auto & first_occurences(ctx.first_occurences);

    std::size_t top{0};
    std::size_t bottom{last_column.size() - 1};

    auto pattern_ix{std::distance(begin, end)};
    assert(pattern_ix >= 0);

    while (top <= bottom)
    {
        if (pattern_ix != 0)
        {
            const auto symbol = begin[-1 + pattern_ix--];
            if (std::find(last_column.begin() + top, last_column.begin() + bottom + 1, symbol) != last_column.begin() + bottom + 1)
            {
                const auto base_ix = ix_by_base(symbol);

                top = first_occurences.at(base_ix) + count.value(base_ix, top, last_column);
                bottom = first_occurences.at(base_ix) + count.value(base_ix, bottom + 1, last_column) - 1;
            }
            else
            {
                return TOP_BOTTOM_INVALID;
            }
        }
        else
        {
            break;
        }
    }

    return {top, bottom};
}


std::vector<pos_type>
better_match(const std::string & pattern, const Context & ctx)
{
    const auto & last_column(ctx.last_column);
    const auto & suffix_array(ctx.partial_suffix_array);
    const auto & count(ctx.count);
    const auto & first_occurences(ctx.first_occurences);

    std::size_t top{0};
    std::size_t bottom{last_column.size() - 1};

    std::size_t pattern_ix{pattern.size()};

    while (top <= bottom)
    {
        if (pattern_ix != 0)
        {
            const auto symbol = pattern[-1 + pattern_ix--];
            if (std::find(last_column.begin() + top, last_column.begin() + bottom + 1, symbol) != last_column.begin() + bottom + 1)
            {
                const auto base_ix = ix_by_base(symbol);

                top = first_occurences.at(base_ix) + count.value(base_ix, top, last_column);
                bottom = first_occurences.at(base_ix) + count.value(base_ix, bottom + 1, last_column) - 1;
            }
            else
            {
                return {};
            }
        }
        else
        {
            break;
        }
    }

    std::vector<pos_type> range;

    for (auto ix = top; ix <= bottom; ++ix)
    {
        range.push_back(suffix_array.value(ix, last_column, count, first_occurences));
    }

    return range;
}


struct CompressedText
{
    enum
    {
        WORD_SZ = 64,
        BITS_PER_SYM = 3,
        SYM_PER_WORD = WORD_SZ / BITS_PER_SYM,
        MASK = (1 << BITS_PER_SYM) - 1,
    };

    inline
    char operator[](std::size_t pos) const
    {
        return (m_bitset[pos / SYM_PER_WORD] >> ((pos % SYM_PER_WORD) * BITS_PER_SYM)) & MASK;
    }

    std::string decompress(std::size_t begin, std::size_t end) const
    {
        std::string result(end - begin, 0);
        std::size_t six = 0;

        for (std::size_t ix{begin}; ix < end; ++ix, ++six)
        {
            result[six] = base_by_ix((*this)[ix]);
        }

        return result;
    }

    //const
    std::size_t m_sz;
    //const
    std::vector<uint64_t> m_bitset;
};


CompressedText compress_text(const std::string & text)
{

    std::vector<uint64_t> bitset((text.size() + CompressedText::SYM_PER_WORD - 1) / CompressedText::SYM_PER_WORD);
    std::size_t bitoff = 0;
    std::size_t word_ix = 0;

    for (const char base : text)
    {
        const auto base_ix = ix_by_base(base);
        bitset[word_ix] |= (base_ix << bitoff);
        if (bitoff == CompressedText::BITS_PER_SYM * (CompressedText::SYM_PER_WORD - 1))
        {
            bitoff = 0;
            ++word_ix;
        }
        else
        {
            bitoff += CompressedText::BITS_PER_SYM;
        }
    }

    return CompressedText{text.size(), bitset};
}


std::vector<std::pair<pos_type, std::size_t>>
approximate_better_match(
    const std::string & pattern,
    const Context & ctx,
    const CompressedText & comptext,
    std::size_t dist)
{
    struct DsetGen
    {
        DsetGen(const std::string & word, const std::size_t dist) :
            m_word(word)
        {
            const auto L = word.size();

            // linspace
            {
                const float step = L / (dist + 1.f);
                m_splits.push_back(0);
                for (std::size_t ix = 1; ix <= dist; ++ix)
                {
                    m_splits.push_back(ix * step);
                }
                m_splits.push_back(L);
            }
        }

        std::pair<std::string, std::size_t> yield()
        {
            if (m_splits.size() >= 2)
            {
                const auto back = m_splits.back();
                m_splits.pop_back();

                return {m_word.substr(m_splits.back(), back - m_splits.back()), m_splits.back()};
            }
            else
            {
                return {"", std::string::npos};
            }
        }

        const std::string & m_word;
        std::vector<std::size_t> m_splits;
    };


    std::set<std::pair<pos_type, std::size_t>> grand_matched;
    auto dset_gen = DsetGen(pattern, dist);
    for (auto dset = dset_gen.yield(); dset.second != std::string::npos; dset = dset_gen.yield())
    {
        for (const auto ix : better_match(dset.first, ctx))
        {
            const auto begin = ix - dset.second;
            const auto end = begin + pattern.size();
            if (begin >= 0 && end <= ctx.last_column.size())
            {
                const auto text = comptext.decompress(begin, end);
                const auto real_dist = std::inner_product(text.cbegin(), text.cend(), pattern.cbegin(),
                    0u,
                    [](const std::size_t acc, const std::size_t val){return acc + val;},
                    [](const char lhs, const char rhs){return lhs != rhs;}
                );
                if (real_dist <= dist)
                {
                    grand_matched.emplace(begin, real_dist);
                }
            }
        }
    }
    // TODO: unordered map

    const std::vector<std::pair<pos_type, std::size_t>> range(grand_matched.cbegin(), grand_matched.cend());
    return range;
//    return grand_matched;
}


typedef uint32_t hash_type;


hash_type kmer_hash_fwd(const char * what, std::size_t depth)
{
    hash_type result{0};

    while (depth--)
    {
        result = (result << 2) | BW::ix_by_nucleobase(*what--);
    }

    return result;
}


std::vector<pos_type>
cached_better_match(
    const char * begin,
    const char * end,
    const Context & ctx,
    const std::unordered_map<hash_type, std::pair<BW::pos_type, BW::pos_type>> & cache,
    std::size_t cache_depth)
{
    const auto hash = kmer_hash_fwd(end - 1, cache_depth);
    if (cache.count(hash) == 0)
    {
        return {};
    }

    const auto & last_column(ctx.last_column);
    const auto & suffix_array(ctx.partial_suffix_array);
    const auto & count(ctx.count);
    const auto & first_occurences(ctx.first_occurences);

    const auto cached_top_bottom = cache.at(hash);

    const auto bw_top_bottom = top_bottom(begin, end - cache_depth, ctx, cached_top_bottom.first, cached_top_bottom.second);

    if (bw_top_bottom == TOP_BOTTOM_INVALID)
    {
        return {};
    }

    std::vector<pos_type> range;

    for (auto ix = bw_top_bottom.first; ix <= bw_top_bottom.second; ++ix)
    {
        range.push_back(suffix_array.value(ix, last_column, count, first_occurences));
    }

    return range;
}


std::vector<std::pair<pos_type, std::size_t>>
cached_approximate_better_match(
    const std::string & pattern,
    const Context & ctx,
    const CompressedText & comptext,
    std::size_t dist,
    const std::unordered_map<hash_type, std::pair<BW::pos_type, BW::pos_type>> & cache,
    std::size_t cache_depth)
{
    struct DsetGen
    {
        DsetGen(const std::string & word, const std::size_t dist) :
            m_word(word)
        {
            const auto L = word.size();

            // linspace
            {
                const float step = L / (dist + 1.f);
                m_splits.push_back(0);
                for (std::size_t ix = 1; ix <= dist; ++ix)
                {
                    m_splits.push_back(ix * step);
                }
                m_splits.push_back(L);
            }
        }

        std::pair<std::string, std::size_t> yield()
        {
            if (m_splits.size() >= 2)
            {
                const auto back = m_splits.back();
                m_splits.pop_back();

                return {m_word.substr(m_splits.back(), back - m_splits.back()), m_splits.back()};
            }
            else
            {
                return {"", std::string::npos};
            }
        }

        const std::string & m_word;
        std::vector<std::size_t> m_splits;
    };

    ////////////////////////////////////////////////////////////////////////////

    std::set<std::pair<pos_type, std::size_t>> grand_matched;

    auto dset_gen = DsetGen(pattern, dist);

    for (auto dset = dset_gen.yield(); dset.second != std::string::npos; dset = dset_gen.yield())
    {
        for (const auto ix : cached_better_match(
            dset.first.c_str(), dset.first.c_str() + dset.first.size(),
            ctx, cache, cache_depth))
        {
            const auto begin = ix - dset.second;
            const auto end = begin + pattern.size();
            if (begin >= 0 && end <= ctx.last_column.size())
            {
                const auto text = comptext.decompress(begin, end);
                const auto real_dist = std::inner_product(text.cbegin(), text.cend(), pattern.cbegin(),
                    0u,
                    [](const std::size_t acc, const std::size_t val){return acc + val;},
                    [](const char lhs, const char rhs){return lhs != rhs;}
                );
                if (real_dist <= dist)
                {
                    grand_matched.emplace(begin, real_dist);
                }
            }
        }
    }
    // TODO: unordered map

    const std::vector<std::pair<pos_type, std::size_t>> range(grand_matched.cbegin(), grand_matched.cend());
    return range;
//    return grand_matched;
}


}  // namespace BW


#endif /* BURROWSWHEELER_HPP_ */
