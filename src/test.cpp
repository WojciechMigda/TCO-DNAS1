/*******************************************************************************
 * Copyright (c) 2016 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: test.cpp
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
 * 2016-04-19   wm              Initial version
 *
 ******************************************************************************/

#include "BW/BurrowsWheeler.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <iostream>
#include <cassert>
#include <string>
#include <algorithm>

using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::StrEq;


TEST(FullSuffixArray, IsConstructedFromPattern_303_4)
{
    const std::string pattern("AATCGGGTTCAATCGGGGT$");
    const std::vector<BW::pos_type> res = BW::full_suffix_array(pattern);
    ASSERT_THAT(res, ElementsAreArray({19, 10, 0, 11, 1, 9, 13, 3, 14, 15, 4, 16, 5, 17, 6, 18, 8, 12, 2, 7}));
}


TEST(FullSuffixArray, IsConstructedFromPattern_304_6)
{
    const std::vector<BW::pos_type> res = BW::full_suffix_array(std::string("ACATGCTACTTT$"));
    ASSERT_THAT(res, ElementsAreArray({12, 0, 7, 2, 1, 5, 8, 4, 11, 6, 3, 10, 9}));
}


TEST(FullSuffixArray, IsConstructedFromPattern_310_2)
{
    const std::vector<BW::pos_type> res = BW::full_suffix_array(std::string("AACGATAGCGGTAGA$"));
    ASSERT_THAT(res, ElementsAreArray({15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5}));
}


TEST(LastColumn, IsConstructedFromPattern303_4)
{
    const std::string pattern("AATCGGGTTCAATCGGGGT$");
    const std::vector<BW::pos_type> full_suffix_array = BW::full_suffix_array(pattern);
    const std::string res = BW::last_column(pattern, full_suffix_array);

    ASSERT_THAT(res, StrEq("TC$AATTTCGCGGGGGTAAG"));
}


TEST(LastColumn, IsConstructedFromPattern304_6)
{
    const std::string pattern("ACATGCTACTTT$");
    const std::vector<BW::pos_type> full_suffix_array = BW::full_suffix_array(pattern);
    const std::string res = BW::last_column(pattern, full_suffix_array);

    ASSERT_THAT(res, StrEq("T$TCAGATTCATC"));
}


int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
