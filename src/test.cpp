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

using ::testing::ContainerEq;
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


TEST(PartialSuffixArray, IsExpandedToFull303_4)
{
    enum {SKIP = 3};

    const std::string pattern("AATCGGGTTCAATCGGGGT$");
    const std::vector<BW::pos_type> full_suffix_array = BW::full_suffix_array(pattern);
    const std::string last_column = BW::last_column(pattern, full_suffix_array);
    const auto count = BW::count(last_column, SKIP);
    const auto first_occurrences = BW::first_occurences(count, last_column);
    const auto partial_suffix_array = BW::partial_suffix_array(full_suffix_array, SKIP);

    std::vector<BW::pos_type> res;
    for (std::size_t ix{0}; ix < full_suffix_array.size(); ++ix)
    {
        res.push_back(partial_suffix_array.value(ix, last_column, count, first_occurrences));
    }
    ASSERT_THAT(res, ContainerEq(full_suffix_array));
}


TEST(CompressedText, IsDecompressed)
{
    const std::string instring = "ACATGCTACTTTAACGATAGCGGTAGAAACGATAGCGGTAGA$";

    const auto compressed = BW::compress_text(instring);

    const auto decompressed = compressed.decompress(4, instring.size() - 4);
    ASSERT_EQ(decompressed, instring.substr(4, instring.size() - 4 - 4));
}


TEST(CompressedText, IsDecompressedFromBeginning)
{
    const std::string instring = "ACATGCTACTTTAACGATAGCGGTAGAAACGATAGCGGTAGA$";

    const auto compressed = BW::compress_text(instring);

    const auto decompressed = compressed.decompress(0, instring.size() - 3);
    ASSERT_EQ(decompressed, instring.substr(0, instring.size() - 3));
}


TEST(CompressedText, IsDecompressedTillEnd)
{
    const std::string instring = "ACATGCTACTTTAACGATAGCGGTAGAAACGATAGCGGTAGA$";

    const auto compressed = BW::compress_text(instring);

    const auto decompressed = compressed.decompress(3, instring.size());
    ASSERT_EQ(decompressed, instring.substr(3, instring.size() - 3));
}


//TEST(ApproximateBWMatch, FindsApproximateMatches)
//{
//    enum {SKIP = 3};
//
//    const std::string pattern("ACATGCTACTTT$");
//    const std::vector<BW::pos_type> full_suffix_array = BW::full_suffix_array(pattern);
//
//    const std::string last_column = BW::last_column(pattern, full_suffix_array);
//    const auto count = BW::count(last_column, SKIP);
//    const auto first_occurrences = BW::first_occurences(count, last_column);
//    const auto partial_suffix_array = BW::partial_suffix_array(full_suffix_array, SKIP);
//    const auto compressed = BW::compress_text(pattern);
//
//    const auto matches = BW::approximate_better_match("ATT",
//        BW::Context{last_column, count, first_occurrences, partial_suffix_array},
//        compressed,
//        1);
//    ASSERT_THAT(matches, ElementsAreArray({2, 7, 8, 9}));
//}


int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
