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

namespace BW
{


typedef uint32_t pos_type;


template <typename SeqT>
std::vector<pos_type>
full_suffix_array(const SeqT & text)
{
    std::vector<pos_type> sarray(text.size());
//    std::vector<index_type> sarray;
//    sarray.resize(text.size());

    std::iota(sarray.begin(), sarray.end(), 0);

    std::sort(sarray.begin(), sarray.end(),
        [&text](const pos_type lhix, const pos_type rhix) -> bool
        {
            typename SeqT::const_pointer lhs = text.data() + lhix;
            typename SeqT::const_pointer rhs = text.data() + rhix;

            const auto result = std::memcmp(lhs, rhs, text.size() - std::max(lhix, rhix));
            return result < 0;
        });

    return sarray;
}


template <typename SeqT>
SeqT last_column(const SeqT & text, const std::vector<uint32_t> & full_suffix_array)
{
//    SeqT last_col(text.size(), {});
    SeqT last_col;
    last_col.resize(text.size());

    assert(full_suffix_array.size() == last_col.size());

    std::transform(full_suffix_array.cbegin(), full_suffix_array.cend(),
        last_col.begin(),
        [&text](const uint32_t ix)
        {
            return ix == 0 ? '$' : text[ix - 1];
        });

    return last_col;
}


constexpr char base_by_ix(const std::size_t ix)
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


constexpr char nucleobase_by_ix(const std::size_t ix)
{
    return ((
        ((uint32_t)'T' << 24) |
        ((uint32_t)'G' << 16) |
        ((uint32_t)'C' << 8) |
        ((uint32_t)'A' << 0) |
        0
        ) >> (ix * 8)) & 0xFF;
}


constexpr std::size_t ix_by_nucleobase(const char base)
{
    return
        ((((0ULL << ('A' - 'A')) |
            (1ULL << ('C' - 'A')) |
            (2ULL << ('G' - 'A')) |
            (3ULL << ('T' - 'A'))) >> (base - 'A')) & 0x3);
}


constexpr std::size_t ix_by_base(const char base)
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
    const std::size_t m_skip;
    const std::vector<std::vector<count_type>> m_counts;

    template <typename SeqT>
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
    std::vector<std::vector<count_type>> result(5);
    std::vector<count_type> stat(5);

    for (auto & vec : result)
    {
        vec.push_back(0);
    }

    for (std::size_t ix{0}; ix < last_column.size(); ++ix)
    {
        const auto item = last_column[ix];

        for (int ix{0}; ix < 5; ++ix)
        {
            stat[ix] += (item == base_by_ix(ix));
        }

        if (ix % SKIP == (SKIP - 1))
        {
            for (int ix{0}; ix < 5; ++ix)
            {
                result[ix].push_back(stat[ix]);
            }
        }
    }

    return Count{SKIP, result};
}


template <typename SeqT>
std::vector<std::size_t>
first_occurences(const Count & count, const SeqT & last_column)
{
    std::vector<std::size_t> first_occ(5);

    first_occ[ix_by_base('$')] = 0;
    first_occ[ix_by_base('A')] = 1;
    first_occ[ix_by_base('C')] = count.value(ix_by_base('A'), last_column.size(), last_column) + first_occ[ix_by_base('A')];
    first_occ[ix_by_base('G')] = count.value(ix_by_base('C'), last_column.size(), last_column) + first_occ[ix_by_base('C')];
    first_occ[ix_by_base('T')] = count.value(ix_by_base('G'), last_column.size(), last_column) + first_occ[ix_by_base('G')];

    return first_occ;
}


struct PartialSuffixArray
{
    const std::size_t m_skip;
    const std::unordered_map<uint32_t, uint32_t> m_sufarr;

    template <typename SeqT>
    uint32_t value(
        uint32_t ix,
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
partial_suffix_array(const std::vector<uint32_t> & full_suffix_array, const std::size_t SKIP = 4)
{
    if (full_suffix_array.size() == 0)
    {
        return PartialSuffixArray{SKIP, {}};
    }

    std::unordered_map<uint32_t, uint32_t> partial;
    for (std::size_t ix{0}; ix < full_suffix_array.size(); ix += SKIP)
    {
        partial[ix] = full_suffix_array[ix];
    }

    return PartialSuffixArray{SKIP, partial};
}


struct Context
{
    const std::string last_column;
    const BW::Count count;
    const std::vector<std::size_t> first_occurences;
    const BW::PartialSuffixArray partial_suffix_array;
};


std::vector<uint32_t>
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

    std::vector<uint32_t> range;

    for (auto ix = top; ix <= bottom; ++ix)
    {
        range.push_back(suffix_array.value(ix, last_column, count, first_occurences));
    }

    return range;
}


}  // namespace BW


#endif /* BURROWSWHEELER_HPP_ */
