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

#include <string>
#include <vector>

struct DNASequencing
{
    int initTest(int testDifficulty)
    {
        return 0;
    }

    int preProcessing()
    {
        return 0;
    }

    int passReferenceGenome(
        int chromatidSequenceId,
        const std::vector<std::string> & chromatidSequence)
    {
        return 0;
    }

    std::vector<std::string> getAlignment(
        int N,
        double normA,
        double normS,
        const std::vector<std::string> & readName,
        const std::vector<std::string> & readSequence)
    {
        std::vector<std::string> ret(N, "");
        for (int i = 0; i < N; ++i)
        {
            std::string qname = "sim" + std::to_string(1 + i / 2) + '/'
                + ((i % 2) ? '2' : '1');
            ret[i] = qname + ",20,1,150,+,0";
        }
        return ret;
    }
};

#endif /* DNA_SEQUENCING_HPP_ */
