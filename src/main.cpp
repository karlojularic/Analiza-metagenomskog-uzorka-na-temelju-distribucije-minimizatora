#include <iostream>
#include <getopt.h>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include "loading_sequences.hpp"
#include "minimizers.hpp"

using namespace std;


auto distribution(const vector<tuple<unsigned int, unsigned int, bool>> minimizers) {
    unordered_map<unsigned int, unsigned int> distribution_vector;

    for (const auto& [kmer, pos, is_original] : minimizers) {
        distribution_vector[kmer]++;
    }

    return distribution_vector;

}


int main(int argc, char *argv[]){
    cout << "Starting program..." << endl;
    unsigned int k = 10;
    unsigned int w = 3;
    
    int opt;
    while ((opt = getopt(argc, argv, "k:w:")) != -1){
        switch(opt) {
            case 'k': k = stoi(optarg); break;
            case 'w': w = stoi(optarg); break;

            default: return 1;
        }
    }

    if (optind + 2 != argc) {
        cerr << "Error: reference and fragments files are expected.\n";
        return 1;
    }

    string file1 = argv[optind];
    string file2 = argv[optind + 1];

    vector<analysis::Sequence> references = analysis::LoadSequences(file1);
    vector<analysis::Sequence> fragments = analysis::LoadSequences(file2);

    for (const auto& ref : references) {
        vector<tuple<unsigned int, unsigned int, bool>> ref_minimizers = analysis::Minimize(ref.seq.c_str(), ref.seq.size(), k, w); 
        unordered_map<unsigned int, unsigned int> distribution_vector = distribution(ref_minimizers);

        cout << "Distribution for reference: " << ref.name << endl;

        for (const auto& [kmer, count] : distribution_vector) {
            cout << kmer << ": " << count << endl;
        }

    }

    
    
    return 0;
}