/*******************************************************************************
 * Copyright (c) 2016 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: main.cpp
 *
 * Description:
 *      TCO-DNAS1
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2016-04-15   wm              Initial version
 *
 ******************************************************************************/

#if 0

>sim11/1
GTCAAATGGTATTTCTGGTTCTACATCCTTGAGGAATCGCCACACTGTCTTCCACAATGGTTGAACTAATTTACATTCCCACCAACAGTGTAAAAGTATTCTTATTTCTCCACAGCCTTGCCAGCATCTATTGTTTCTTGACTTTTTAAT

>sim11/2
ATAAAAACCCTAGAAGAAAACCTAGGCAATACCATTTAGGACATAGGCATGGGCAAAGACTTCATGAAGAAAATACCAAAAGCAATTGCAACAAAAGCCAAAATTGTCAAATGGGATCTAATTAAATGAACGAGCTTCTGCACAGCAAAA

TTTTGCTGTGCAGAAGCTCGTTCATTTAATTAGATCCCATTTGACAATTTTGGCTTTTGTTGCAATTGCTTTTGGTATTTTCTTCATGAAGTCTTTGCCCATGCCTATGTCCTAAATGGTATTGCCTAGGTTTTCTTCTAGGGTTTTTAT

sim11/1,20,98162,98311,+,GTCAAATGGTATTTCTGGTTCTACATCCTTGAGGAATCGCCACACTGTCTTCCACAATGGTTGAACTAATTTACATTCCCACCAACAGTGTAAAAGTATTCTTATTTCTCCACAGCCTTGCCAGCATCTATTGTTTCTTGACTTTTTAAT
sim11/2,20,98594,98743,-,TTTTGCTGTGCAGAAGCTCGTTCATTTAATTAGATCCCATTTGACAATTTTGGCTTTTGTTGCAATTGCTTTTGGTATTTTCTTCATGAAGTCTTTGCCCATGCCTATGTCCTAAATGGTATTGCCTAGGTTTTCTTCTAGGGTTTTTAT


#endif

#include "dna_sequencing.hpp"

#include <boost/program_options.hpp>
#include <boost/any.hpp>

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <exception>
#include <string>
#include <tuple>
#include <fstream>
#include <cassert>

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::variables_map args;

    {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("test-type", po::value<int>()->default_value(0), "test type, 0: small, 1: medium, 2: large")
            ("test-number", po::value<int>()->default_value(5), "test number, [5-10] for small and medium, 5 for large")
            ;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc), args);
        }
        catch (po::error & ex)
        {
            std::cout << ex.what() << std::endl;
        }
        catch (boost::bad_any_cast & ex)
        {
            std::cout << ex.what() << std::endl;
        }
        po::notify(args);

        if (args.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }
    }

    const int test_type = std::max(std::min(args["test-type"].as<int>(), 2), 0);
    const int test_number = test_type < 2 ? std::max(std::min(args["test-number"].as<int>(), 10), 5) : 5;

    std::cout << "[main] Executing test of difficulty " << test_type << std::endl;
    std::cout << "[main]           test number: " << test_number << std::endl;

    const std::tuple<std::string, std::vector<int>, double, double> test_config[] =
    {
        std::make_tuple("small", std::vector<int>{20},
                        -3.392, 0.5),
        std::make_tuple("medium", std::vector<int>{1, 11, 20},
                        -3.962, 0.5),
        std::make_tuple("large", std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                        -2.710, 0.5),
    };

    const std::string input_prefix = "../data/" + std::get<0>(test_config[test_type]) + std::to_string(test_number) + ".";
    const std::string chroma_prefix = "../data/chromatid";

    {
        DNASequencing dna_sequencing;

        dna_sequencing.initTest(test_type);


        for (const int chroma_id : std::get<1>(test_config[test_type]))
        {
            std::cout << "[main] reading chromatid " << chroma_id << std::endl;

            const std::string chroma_path = chroma_prefix + std::to_string(chroma_id) + ".fa";
            std::ifstream ifs(chroma_path);
            std::string line;

            std::getline(ifs, line); // skip header

            std::vector<std::string> chromatid;
            while (std::getline(ifs, line))
            {
                if (line.back() == '\r')
                {
                    line.pop_back();
                }
                chromatid.push_back(line);
            }
            dna_sequencing.passReferenceGenome(chroma_id, chromatid);
        }

        dna_sequencing.preProcessing();

        const std::string fa1_path = input_prefix + "fa1";
        const std::string fa2_path = input_prefix + "fa2";

        std::cout << "[main] reading read pairs" << std::endl;

        std::ifstream ifs_fa1(fa1_path);
        std::ifstream ifs_fa2(fa2_path);
        std::string line_fa1;
        std::string line_fa2;
        std::vector<std::string> read_pairs;
        std::vector<std::string> read_pair_ids;
        while (std::getline(ifs_fa1, line_fa1) && std::getline(ifs_fa2, line_fa2))
        {
            if (line_fa1.back() == '\r')
            {
                line_fa1.pop_back();
            }
            if (line_fa2.back() == '\r')
            {
                line_fa2.pop_back();
            }
            read_pair_ids.push_back(line_fa1.substr(1, line_fa1.size() - 1));
            read_pair_ids.push_back(line_fa2.substr(1, line_fa2.size() - 1));

            if (std::getline(ifs_fa1, line_fa1) && std::getline(ifs_fa2, line_fa2))
            {
                if (line_fa1.back() == '\r')
                {
                    line_fa1.pop_back();
                }
                if (line_fa2.back() == '\r')
                {
                    line_fa2.pop_back();
                }
                read_pairs.push_back(line_fa1);
                read_pairs.push_back(line_fa2);
            }
            else
            {
                read_pair_ids.pop_back();
                read_pair_ids.pop_back();
            }
        }
        assert(read_pairs.size() == read_pair_ids.size());
        std::cout << "[main] read " << read_pairs.size() << " read pairs from " << input_prefix << "fa{1,2}" << std::endl;

        std::cout << "[main] requesting alignment" << std::endl;

        const std::vector<std::string> results =
            dna_sequencing.getAlignment(read_pairs.size(),
                std::get<2>(test_config[test_type]), std::get<3>(test_config[test_type]),
                read_pair_ids, read_pairs);
    }

    return 0;
}
