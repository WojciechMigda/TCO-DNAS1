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
#include <boost/algorithm/string/split.hpp>

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <exception>
#include <string>
#include <tuple>
#include <fstream>
#include <cassert>
#include <utility>
#include <vector>
#include <map>


boost::program_options::variables_map
parse_options(int argc, char **argv)
{
    namespace po = boost::program_options;

    po::variables_map args;

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
        std::exit(1);
    }

    return args;
}


std::vector<std::string>
chromatid_from(const std::string & chroma_path)
{
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

    return chromatid;
}


std::pair<std::vector<std::string>, std::vector<std::string>>
read_pairs_from(
    const std::string & fa1_path,
    const std::string & fa2_path)
{
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

    return {read_pair_ids, read_pairs};
}


double score(
    const std::vector<std::string> & results,
    const double norm_a,
    const std::string & minisam_path)
{
    struct Position
    {
        int rname;
        int from;
        int to;
        char strand;
    };

    struct ReadResult
    {
        double confidence;
        int r;
    };

    auto parse_truth = [](const std::string & path)
    {
        std::map<std::string, Position> res;
        std::ifstream ifs(path);
        std::string s;
        while (ifs >> s)
        {
            std::vector<std::string> ovec;
            boost::split(ovec, s, [](char c){return c == ',';}, boost::token_compress_off);
            res[ovec[0]] = Position{std::stoi(ovec[1]), std::stoi(ovec[2]), std::stoi(ovec[3]), ovec[4][0]};
        }

        return res;
    };

    std::map<std::string, Position> truth = parse_truth(minisam_path);

    constexpr int MAX_POSITION_DIST = 300;
    auto build_read_results = [](const std::map<std::string, Position>& truth,
        const std::vector<std::string>& results)
    {
        std::vector<ReadResult> read_results;
        const int n = results.size();
        int correct = 0;
        for (int i = 0; i < n; ++i)
        {
            std::vector<std::string> tokens;
            boost::split(tokens, results[i], [](char c){return c == ',';}, boost::token_compress_off);
            const auto p = truth.find(tokens[0]);
            const Position & position = p->second;
            int r = 1;
            r = (std::stoi(tokens[1]) == position.rname) ? r : 0;
            r = (tokens[4][0] == position.strand) ? r : 0;
            const int start0 = stoi(tokens[2]);
            const int start1 = position.from;
            r = (std::abs(start0 - start1) < MAX_POSITION_DIST) ? r : 0;
            const double confidence = std::stod(tokens[5]);
            read_results.push_back(ReadResult{confidence, r});
            correct += r;
        }
        std::cout << "Number of correct answers: " << correct << '/' << n << " = "
            << (double)correct / (double)n << std::endl;

        return read_results;
    };

    std::vector<ReadResult> read_results = build_read_results(truth, results);

    auto compute_accuracy = [](std::vector<ReadResult> & read_results, double norm_a)
    {
        const int n = read_results.size();
        std::sort(read_results.begin(), read_results.end(),
            [](const ReadResult& lhs, const ReadResult& rhs)
            {
                return lhs.confidence > rhs.confidence;
            });
        // merge results of equal confidence
        std::vector<int> cumul_si{read_results[0].r};
        std::vector<int> pos{0};
        for (int i = 1; i < n; ++i)
        {
            if (read_results[i].confidence == read_results[i - 1].confidence)
            {
                cumul_si.back() += read_results[i].r;
                pos.back() = i;
            }
            else
            {
                const double cumul = cumul_si.back() + read_results[i].r;
                cumul_si.push_back(cumul);
                pos.push_back(i);
            }
        }
        // compute the AuC
        double auc = 0.0;
        const double invn = 1.0 / (double)n;
        const double invnp1 = 1.0 / (double)(n + 1);
        const double lfmultiplier = 1.0 / std::log(n + 1);
        const int m = cumul_si.size();
        for (int i = 0; i < m; ++i)
        {
            const double fi = 1.0 * (2 + pos[i] - cumul_si[i]) * invnp1;
            const double fi1 =
                (i == (m - 1)) ?
                    1.0 : 1.0 * (2 + pos[i + 1] - cumul_si[i + 1]) * invnp1;
            const double lfi = lfmultiplier * log(fi);
            const double lfi1 = lfmultiplier * log(fi1);
            auc += cumul_si[i] * (lfi1 - lfi) * invn;
        }
        std::cout << "auc = " << auc << std::endl;

        constexpr double MAX_AUC = 0.999999;
        const double tmp = std::log(1 - std::min(auc, MAX_AUC));
        std::cout << "log(1 - min(auc, MAX_AUC)) = " << tmp << std::endl;
        std::cout << "NormA = " << norm_a << std::endl;
        const double accuracy = tmp / norm_a;
        std::cout << "accuracy = " << accuracy << std::endl;

        return accuracy;
    };

    return compute_accuracy(read_results, norm_a);
}


int main(int argc, char **argv)
{
    const auto args = parse_options(argc, argv);


    const int test_type = std::max(std::min(args.at("test-type").as<int>(), 2), 0);
    const int test_number = test_type < 2 ? std::max(std::min(args.at("test-number").as<int>(), 10), 5) : 5;

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

            const std::vector<std::string> chromatid = chromatid_from(chroma_path);

            dna_sequencing.passReferenceGenome(chroma_id, chromatid);
        }

        dna_sequencing.preProcessing();

        const std::string fa1_path = input_prefix + "fa1";
        const std::string fa2_path = input_prefix + "fa2";

        std::cout << "[main] reading read pairs" << std::endl;

        const auto read_pairs_w_ids = read_pairs_from(fa1_path, fa2_path);
        const auto & read_pair_ids = read_pairs_w_ids.first;
        const auto & read_pairs = read_pairs_w_ids.second;

        assert(read_pairs.size() == read_pair_ids.size());
        std::cout << "[main] read " << read_pairs.size() << " read pairs from " << input_prefix << "fa{1,2}" << std::endl;

        std::cout << "[main] requesting alignment" << std::endl;

        const std::vector<std::string> results =
            dna_sequencing.getAlignment(read_pairs.size(),
                std::get<2>(test_config[test_type]), std::get<3>(test_config[test_type]),
                read_pair_ids, read_pairs);

        std::cout << results.front() << std::endl;

        const std::string minisam_path = input_prefix + "minisam";
        const double SCORE = score(results, std::get<2>(test_config[test_type]), minisam_path);

        std::cout << "SCORE: " << SCORE << std::endl;
    }

    return 0;
}
