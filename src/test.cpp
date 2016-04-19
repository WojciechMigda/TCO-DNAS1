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

TEST(SuffixArray, IsConstructedFromPattern)
{
    const std::vector<BW::pos_type> res = BW::full_suffix_array(std::string("ACATGCTACTTT$"));
    ASSERT_THAT(res, ElementsAreArray({12, 0, 7, 2, 1, 5, 8, 4, 11, 6, 3, 10, 9}));
}


int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
