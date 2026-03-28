#include <iostream>
#include <getopt.h>
#include "bioparser/fasta_parser.hpp"
#include <vector>
#include <tuple>
#include <string>

using namespace std;

inline unsigned char charTo2Bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

inline unsigned int getReverseComplement(unsigned int kmer, unsigned int k) {
    unsigned int rc = 0;
    for (unsigned int i = 0; i < k; ++i) {
        unsigned int baza = kmer & 0x03; 
    
        unsigned int komplement = 3 - baza;
        
        rc = (rc << 2) | komplement;

        kmer >>= 2;
    }
    return rc;
}

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

    const char* sequence;
    unsigned int sequence_len;
    unsigned int kmer_len;
    unsigned int window_len;

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

    if (sequence_len < kmer_len) return minimizers;

    int num_kmers = sequence_len - kmer_len + 1;

    struct KmerInfo {
        unsigned int value;
        unsigned int pos;
        bool is_original;
    };
    std::vector<KmerInfo> canonical_kmers;

    for (int i = 0; i < num_kmers; ++i) {
        unsigned int current_kmer = 0;
        for (unsigned int j = 0; j < kmer_len; ++j) {
            current_kmer = (current_kmer << 2) | charTo2Bit(sequence[i + j]);
        }
        
        unsigned int rc_kmer = getReverseComplement(current_kmer, kmer_len);

        if (current_kmer < rc_kmer) {
            canonical_kmers.push_back({current_kmer, int(i), true});
        } else {
            canonical_kmers.push_back({rc_kmer, int(i), false});
        }
    }

    if (canonical_kmers.size() < window_len){
        window_len = canonical_kmers.size();
    }

    // za izbjegavanje dodavanja istog minimizera
    int last_added_pos = -1;

    for (int i = 0; i <= (int)(canonical_kmers.size() - window_len); ++i) {
        KmerInfo min_in_window = canonical_kmers[i];

        for(unsigned int j = 1; j < window_len; ++j) {
            if (canonical_kmers[i + j].value < min_in_window.value) {
                min_in_window = canonical_kmers[i + j];
            }
        }

        if (last_added_pos != min_in_window.pos) {
            minimizers.push_back(std::make_tuple(
                min_in_window.value,
                min_in_window.pos,
                min_in_window.is_original
            ));
            last_added_pos = (int)min_in_window.pos;
        }
    }

    return minimizers;

}

int main(int argc, char *argv[]){
    cout << "Starting program..." << endl;
    /* unsigned int k = 10;
    unsigned int w = 3;
    
    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "k:w:", &option_index)) != -1){
        swith(opt) {
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
    string file2 = argv[optind + 1]; */

    auto minimizers = parsingFasta();

    for (const auto &minimizer : minimizers) {
        cout << "Minimizer: " << std::get<0>(minimizer) << ", Position: " << std::get<1>(minimizer) << ", Is Original: " << std::get<2>(minimizer) << endl;
    }

    return 0;
}