#include <iostream>
#include <getopt.h>
#include "bioparser/fasta_parser.hpp"

using namespace std;

auto parsingFasta(){

    struct Sequence {
        string name;
        string seq;

        public:
        Sequence(
            const char* n, uint32_t n_len,
            const char* s, uint32_t s_len) 
            : name(n, n_len), seq(s, s_len) {}
        };

    auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>("../data/medical-waste-metagenome-blaOXA.fasta");
    auto s = p->Parse(-1);

    vector<Sequence> references;
    for (auto &seq_ptr : s){
        references.push_back(*seq_ptr);
    }

    for (auto &seq : references) {
        cerr << seq.name << " : " << seq.seq.length() << "\n";
    }

}

int main(){
    cout << "Starting program..." << endl;

    parsingFasta();

    return 0;
}