#include "minimizers.hpp"

namespace analysis {

std::vector<std::pair<unsigned int, unsigned int>> Minimize(
    const char* sequence, 
    unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len) {

    std::vector<std::pair<unsigned int, unsigned int>> minimizers;

    if (sequence_len < kmer_len) return minimizers;

    int num_kmers = sequence_len - kmer_len + 1;
    
    struct KmerInfo {
        unsigned int value;
        unsigned int pos;
    };
    std::vector<KmerInfo> kmers;

    for (int i = 0; i < num_kmers; ++i) {
        unsigned int current_kmer = 0;
        for (unsigned int j = 0; j < kmer_len; ++j) {
            current_kmer = (current_kmer << 2) | charTo2Bit(sequence[i + j]);
        }

        kmers.push_back({current_kmer, (unsigned int)i});

    }

    if (kmers.size() < window_len) {
        window_len = kmers.size();
    }

    int last_added_pos = -1;

    for (int i = 0; i <= (int)kmers.size() - (int)window_len; ++i) {
        KmerInfo min_in_window = kmers[i];

        for (unsigned int j = 1; j < window_len; ++j) {
            if (kmers[i + j].value < min_in_window.value) {
                min_in_window = kmers[i + j];
            }
        }

        if (last_added_pos != (int)min_in_window.pos) {
            minimizers.push_back(std::make_pair(
                min_in_window.value, 
                min_in_window.pos
            ));
            last_added_pos = (int)min_in_window.pos;
        }
    }

    return minimizers;
}

}