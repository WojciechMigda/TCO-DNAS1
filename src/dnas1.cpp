#include <bits/stdc++.h>
using namespace std;

// If you submit, comment the following line
#define LOCAL_ONLY

class DNASequencing
{
public:
    int initTest(int testDifficulty)
    {
        return 0;
    }

    int preProcessing()
    {
        return 0;
    }

    int passReferenceGenome(int chromatidSequenceId,
        const vector<string>& chromatidSequence)
    {
        return 0;
    }

    vector<string> getAlignment(int N, double normA, double normS,
        const vector<string>& readName, const vector<string>& readSequence)
    {
        vector<string> ret(N, "");
        for (int i = 0; i < N; ++i)
        {
            string qname = "sim" + to_string(1 + i / 2) + '/'
                + ((i % 2) ? '2' : '1');
            ret[i] = qname + ",20,1,150,+,0";
        }
        return ret;
    }
};

#ifdef LOCAL_ONLY
/**
 * Constants from the problem statement
 */
const int MAX_POSITION_DIST = 300;
const double NORM_A_SMALL = -3.392;
const double NORM_A_MEDIUM = -3.962;
const double NORM_A_LARGE = -2.710;
const double MAX_AUC = 0.999999;

/**
 * Position: describe the position of a read within the genome
 */
struct Position
{
    int rname;
    int from;
    int to;
    char strand;
};

/**
 * ReadResult: result of a read alignment
 */
struct ReadResult
{
    double confidence;
    int r;
};

/**
 * Split a comma-separated string into a vector of string
 * @param row	the string to be split
 * @return	the vector of string
 */
vector<string> tokenize(const string& row)
{
    vector<string> tokens;
    for (int i = 0, pos = 0, n = row.size(); i < n; ++i)
    {
        if (i == n - 1 || row[i + 1] == ',')
        {
            string token = row.substr(pos, (i + 1) - pos);
            tokens.push_back(token);
            pos = i + 2;
        }
    }
    return tokens;
}

/**
 * Read a minisam file and build a map of ground truth
 * @param path	the path of the minisam file storing the ground truth
 * @return a map[read_name] = read_Position
 */
map<string, Position> parse_truth(const string& path)
{
    map<string, Position> res;
    ifstream ifs(path);
    string s;
    while (ifs >> s)
    {
        vector<string> tokens = tokenize(s);
        try
        {
            string qname = tokens[0];
            int chromatid = stoi(tokens[1]);
            int from = stoi(tokens[2]);
            int to = stoi(tokens[3]);
            char strand = tokens[4][0];
            res[qname] = Position
                { chromatid, from, to, strand };
        } catch (exception& e)
        {
            ;
        }
    }
    return res;
}

/**
 * For each string of the results vector, build a read result {confidence, r}
 * @param truth		the map of ground truth position for each read
 * @param results	the vector of results as return by getAlignment
 * @return a vector of ReadResult, that is {confidence, r}
 */
vector<ReadResult> build_read_results(const map<string, Position>& truth,
    const vector<string>& results)
{
    vector<ReadResult> read_results;
    int n = results.size();
    int correct = 0;
    for (int i = 0; i < n; ++i)
    {
        vector<string> tokens = tokenize(results[i]);
        auto p = truth.find(tokens[0]);
        const Position& position = p->second;
        int r = 1;
        r = (stoi(tokens[1]) == position.rname) ? r : 0;
        r = (tokens[4][0] == position.strand) ? r : 0;
        int start0 = stoi(tokens[2]);
        int start1 = position.from;
        r = (abs(start0 - start1) < MAX_POSITION_DIST) ? r : 0;
        double confidence = stod(tokens[5]);
        read_results.push_back(ReadResult
            { confidence, r });
        correct += r;
    }
    cerr << "Number of correct answers: " << correct << '/' << n << " = "
        << (double) correct / (double) n << endl;
    return read_results;
}

/**
 * Compute the accuracy given the {confidence, r} pairs and the normalization facto
 * @param read_results	a vector of {confidence, r} results
 * @param norm_a		as described in the problem statement
 * @return	a double, the computed accuracy
 */
double compute_accuracy(vector<ReadResult>& read_results, double norm_a)
{
    int n = read_results.size();
    sort(read_results.begin(), read_results.end(),
        [](const ReadResult& lhs, const ReadResult& rhs)
        {   return (lhs.confidence>rhs.confidence);});
    // merge results of equal confidence
    vector<int> cumul_si
        { read_results[0].r };
    vector<int> pos
        { 0 };
    for (int i = 1; i < n; ++i)
    {
        if (read_results[i].confidence == read_results[i - 1].confidence)
        {
            cumul_si.back() += read_results[i].r;
            pos.back() = i;
        }
        else
        {
            double cumul = cumul_si.back() + read_results[i].r;
            cumul_si.push_back(cumul);
            pos.push_back(i);
        }
    }
    // compute the AuC
    double auc = 0.0;
    double invn = 1.0 / (double) n;
    double invnp1 = 1.0 / (double) (n + 1);
    double lfmultiplier = 1.0 / log(n + 1);
    int m = cumul_si.size();
    for (int i = 0; i < m; ++i)
    {
        double fi = 1.0 * (2 + pos[i] - cumul_si[i]) * invnp1;
        double fi1 =
            (i == m - 1) ?
                1.0 : 1.0 * (2 + pos[i + 1] - cumul_si[i + 1]) * invnp1;
        double lfi = lfmultiplier * log(fi);
        double lfi1 = lfmultiplier * log(fi1);
        auc += cumul_si[i] * (lfi1 - lfi) * invn;
    }
    cout << "auc = " << auc << endl;
    double tmp = log(1 - min(auc, MAX_AUC));
    cout << "log(1 - min(auc, MAX_AUC)) = " << tmp << endl;
    cout << "NormA = " << norm_a << endl;
    double accuracy = tmp / norm_a;
    cout << "accuracy = " << accuracy << endl;
    return accuracy;
}

/**
 * Perform a single test
 * @param testDifficulty	define the test type (SMALL=0, MEDIUM=1, LARGE=2)
 * @return	alignments in format specified in the problem statement
 */
vector<string> perform_test(int testDifficulty, double norm_a)
{
    // test data path and description
    string fa1_path, fa2_path;
    vector<int> chr_ids;
    if (testDifficulty == 0)
    {
        fa1_path = "../data/small5.fa1";
        fa2_path = "../data/small5.fa2";
        chr_ids = vector<int>
            { 20 };
    }
    else if (testDifficulty == 1)
    {
        fa1_path = "../data/medium5.fa1";
        fa2_path = "../data/medium5.fa2";
        chr_ids = vector<int>
            { 1, 11, 20 };
    }
    else if (testDifficulty == 2)
    {
        fa1_path = "../data/large5.fa1";
        fa2_path = "../data/large5.fa2";
        for (int i = 1; i <= 24; ++i)
            chr_ids.push_back(i);
    }
    // call the MM DNASequencing methods
    DNASequencing dna_sequencing;
    dna_sequencing.initTest(testDifficulty);
    // load chromatid
    for (int chromatid_seq_id : chr_ids)
    {
        vector<string> chromatid_seq;
        string path = "../data/chromatid" + to_string(chromatid_seq_id) + ".fa";
        ifstream ifs(path);
        string s;
        // skip header
        getline(ifs, s);
        cerr << "Skip header: " << s << endl;
        // pack all lines in chromatid_seq
        for (int i = 0; getline(ifs, s); ++i)
        {
            if (s.back() == '\r')
                s.pop_back();
            chromatid_seq.push_back(s);
        }
        dna_sequencing.passReferenceGenome(chromatid_seq_id, chromatid_seq);
    }
    dna_sequencing.preProcessing();
    // load reads
    vector<string> read_id, read_seq;
    {
        ifstream ifs1(fa1_path);
        ifstream ifs2(fa2_path);
        string s1, s2;
        while (getline(ifs1, s1) && getline(ifs2, s2))
        {
            if (s1.back() == '\r')
                s1.pop_back();
            if (s2.back() == '\r')
                s2.pop_back();
            read_id.push_back(s1.substr(1, s1.size() - 1));
            read_id.push_back(s2.substr(1, s2.size() - 1));
            getline(ifs1, s1);
            getline(ifs2, s2);
            if (s1.back() == '\r')
                s1.pop_back();
            if (s2.back() == '\r')
                s2.pop_back();
            read_seq.push_back(s1);
            read_seq.push_back(s2);
        }
    }
    int nreads = read_id.size();
    // compute alignments
    vector<string> results = dna_sequencing.getAlignment(nreads, norm_a, 0.5,
        read_id, read_seq);
    return results;
}

/**
 * Main function: read the data, perform the DNA alignments and score results
 */
int main()
{
    const int testDifficulty = 1;
    string minisam_path;
    double norm_a;
    if (testDifficulty == 0)
    {
        minisam_path = "../data/small5.minisam";
        norm_a = NORM_A_SMALL;
    }
    else if (testDifficulty == 1)
    {
        minisam_path = "../data/medium5.minisam";
        norm_a = NORM_A_MEDIUM;
    }
    else if (testDifficulty == 2)
    {
        minisam_path = "../data/large5.minisam";
        norm_a = NORM_A_LARGE;
    }
    // perform test
    vector<string> results = perform_test(testDifficulty, norm_a);
    // load truth
    map<string, Position> truth = parse_truth(minisam_path);
    vector<ReadResult> read_results = build_read_results(truth, results);
    // scoring
    double accuracy = compute_accuracy(read_results, norm_a);
    return 0;
}
#endif
