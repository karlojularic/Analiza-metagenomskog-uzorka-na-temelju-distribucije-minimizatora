#include <string>
#include <vector>
#include <cstdint>

namespace analysis {

struct Sequence {
    std::string name;
    std::string seq;
    std::string qual;

    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len)
        : name(n, n_len), seq(s, s_len) {}

    Sequence(const char* n, std::uint32_t n_len,
             const char* s, std::uint32_t s_len,
             const char* q, std::uint32_t q_len)
        : name(n, n_len), seq(s, s_len), qual(q, q_len) {}
};

std::vector<Sequence> LoadSequences(const std::string& path);

}