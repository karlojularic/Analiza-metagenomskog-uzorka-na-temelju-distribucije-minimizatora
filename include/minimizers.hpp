#pragma once
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

namespace analysis {

std::vector<std::pair<unsigned int, unsigned int>> Minimize(
    const char* sequence, 
    unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

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

inline std::string getReverseChain(std::string sequence) {
    for (int i = 0; i < sequence.size(); ++i) {
        if (sequence[i] == 'A' || sequence[i] == 'a') sequence[i] = 'T';
        else if (sequence[i] == 'C' || sequence[i] == 'c') sequence[i] = 'G';
        else if (sequence[i] == 'G' || sequence[i] == 'g') sequence[i] = 'C';
        else if (sequence[i] == 'T' || sequence[i] == 't') sequence[i] = 'A';
    }
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

}