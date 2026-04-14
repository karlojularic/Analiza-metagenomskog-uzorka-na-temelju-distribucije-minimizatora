#include <iostream>
#include <getopt.h>
#include <vector>
#include <tuple>
#include <string>
#include "loading_sequences.hpp"
#include "minimizers.hpp"

using namespace std;

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

    vector<tuple<unsigned int, unsigned int, bool>> minimizers = analysis::Minimize(references[0].seq.c_str(), references[0].seq.size(), k, w);

    return 0;
}