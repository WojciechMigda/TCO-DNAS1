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

namespace BW
{


template <typename SeqT>
std::vector<uint32_t>
full_suffix_array(const SeqT & text)
{
    typedef uint32_t index_type;

    std::vector<index_type> sarray(text.size());
//    std::vector<index_type> sarray;
//    sarray.resize(text.size());

    std::iota(sarray.begin(), sarray.end(), 0);

    std::sort(sarray.begin(), sarray.end(),
        [&text](const index_type lhix, const index_type rhix) -> bool
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
        ((uint64_t)'T' << 32) |
        ((uint64_t)'G' << 24) |
        ((uint64_t)'C' << 16) |
        ((uint64_t)'A' << 8) |
        ((uint64_t)'$' << 0)
        ) >> (ix * 8)) & 0xFF;
}


constexpr std::size_t ix_by_base(const char base)
{
    return
        base == '$' ?
            0
            :
            ((((0ULL << ('A' - 'A')) |
               (1ULL << ('C' - 'A')) |
               (2ULL << ('G' - 'A')) |
               (3ULL << ('T' - 'A'))) >> (base - 'A')) & 0x3) + 1;
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

    first_occ[0] = 0; // $
    first_occ[1] = 1; // A
    first_occ[2] = count.value(1, last_column.size(), last_column) + 1;            // C
    first_occ[3] = count.value(2, last_column.size(), last_column) + first_occ[2]; // G
    first_occ[4] = count.value(3, last_column.size(), last_column) + first_occ[3]; // T

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
        std::size_t backtrack{0};

        while (m_sufarr.count(ix) == 0)
        {
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


}  // namespace BW


#endif /* BURROWSWHEELER_HPP_ */
