#include "loading_sequences.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <vector>

namespace analysis {

    std::vector<Sequence> LoadSequences(const std::string& path) {
    std::vector<Sequence> sequences;
    
    bool is_fastq = (path.find(".fastq") != std::string::npos);

    std::unique_ptr<bioparser::Parser<Sequence>> parser;
    if (is_fastq) {
        parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(path);
    } else {
        parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);
    }

    auto parsed_seqs = parser->Parse(-1);
    for (auto& seq_ptr : parsed_seqs) {
        sequences.push_back(*seq_ptr);
    }
    return sequences;

    }

}