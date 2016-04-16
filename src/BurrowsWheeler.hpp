/*******************************************************************************
 * Copyright (c) 2016 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the GNU LGPL v3
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

namespace BW
{

template <typename SeqT>
std::vector<uint32_t>
full_suffix_array(const SeqT & text)
{
    typedef uint32_t index_type;

    std::vector<index_type> sarray;
    sarray.resize(text.size());

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


template <typename SeqT>
std::vector<std::vector<uint32_t>> count(const SeqT & last_column, const std::size_t SKIP = 100)
{
    std::vector<std::vector<uint32_t>> result(5);
    std::vector<uint32_t> stat(5);

    for (auto & vec : result)
    {
        vec.push_back(0);
    }

    for (std::size_t ix{0}; ix < last_column.size(); ++ix)
    {
        const auto item = last_column[ix];

        stat[0] += (item == '$');
        stat[1] += (item == 'A');
        stat[2] += (item == 'C');
        stat[3] += (item == 'G');
        stat[4] += (item == 'T');

        if (ix % SKIP == (SKIP - 1))
        {
            result[0].push_back(stat[0]);
            result[1].push_back(stat[1]);
            result[2].push_back(stat[2]);
            result[3].push_back(stat[3]);
            result[4].push_back(stat[4]);
        }
    }

//    for (const auto item : last_column)
//    {
//        result[0].push_back(result[0].back() + (item == '$'));
//        result[1].push_back(result[1].back() + (item == 'A'));
//        result[2].push_back(result[2].back() + (item == 'C'));
//        result[3].push_back(result[3].back() + (item == 'G'));
//        result[4].push_back(result[4].back() + (item == 'T'));
//    }

    return result;
}


}  // namespace BW


#endif /* BURROWSWHEELER_HPP_ */
