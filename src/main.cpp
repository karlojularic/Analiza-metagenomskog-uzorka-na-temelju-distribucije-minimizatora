#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <map>
#include "loading_sequences.hpp"
#include "minimizers.hpp"

using namespace std;

auto distribution(const vector<tuple<unsigned int, unsigned int, bool>>& minimizers) {
    unordered_map<unsigned int, unsigned int> distribution_vector;

    for (const auto& [kmer, pos, is_original] : minimizers) {
        distribution_vector[kmer]++;
    }

    return distribution_vector;

}

void create_metagenomic_reference(const vector<string>& reference_files, const string& output_file) {
    ofstream out(output_file);

    for (const auto& file : reference_files) {
        ifstream in("../data/" + file);
        string line;
        while (getline(in, line)) {
            out << line << endl;
        }
        in.close();
    }

    out.close();
}

double cosine_similarity(const unordered_map<unsigned int, unsigned int>& ref, const unordered_map<unsigned int, unsigned int>& frag){
    double dot = 0.0;
    double denom_a = 0.0;
    double denom_b = 0.0;
    for (const auto& [kmer, count] : ref) {
        if (frag.find(kmer) != frag.end()) {
            dot += count * frag.at(kmer);
        }
        denom_a += count * count;
    }

    for (const auto& [kmer, count] : frag) {
        denom_b += count * count;
    }

    return dot / (sqrt(denom_a) * sqrt(denom_b));

}

void export_to_csv(const map<string, unordered_map<unsigned int, unsigned int>>& distributions, const string& output_file, unsigned int k) {
    ofstream out(output_file);

    for (const auto& [name, dist] : distributions) {
        out << name;
        if (&name != &distributions.rbegin()->first) {
            out << ";";
        }
    }
    out << endl;

    for (unsigned int i = 0; i < int(pow(4, k)); ++i) {
        for (const auto& [name, dist] : distributions) {
            unsigned int temp = 0;
            if (dist.find(i) != dist.end()) {
                temp = dist.at(i);
            } else {
                temp = 0;
            }

            out << temp;
            if (&name != &distributions.rbegin()->first) {
                out << ";";
            }
        }
        out << endl;
    }
    out.close();
}   


int main(int argc, char *argv[]){
    cout << "Starting program..." << endl;

    string configuration = "../config.txt";
    string lines;
    unsigned int k = 3;
    unsigned int w = 10;
    vector<string> reference_files;
    vector<string> fragment_files;

    if(!(filesystem::exists(configuration))) {
        cerr << "Configuration file does not exist: " << configuration << endl;
        return 1;
    }

    ifstream config_file(configuration);
    while (getline(config_file, lines)) {
        if (lines.rfind("reference_files=", 0) == 0) {
            string files_str = lines.substr(16);
            size_t pos = 0;
            while ((pos = files_str.find(',')) != string::npos) {
                reference_files.push_back(files_str.substr(0, pos));
                files_str.erase(0, pos + 1);
            }
            reference_files.push_back(files_str);
        } else if (lines.rfind("fragment_files=", 0) == 0) {
            string files_str = lines.substr(15);
            size_t pos = 0;
            while ((pos = files_str.find(',')) != string::npos) {
                fragment_files.push_back(files_str.substr(0, pos));
                files_str.erase(0, pos + 1);
            }
            fragment_files.push_back(files_str);
        } else if (lines.rfind("k=", 0) == 0) {
            k = stoi(lines.substr(2));
        } else if (lines.rfind("w=", 0) == 0) {
            w = stoi(lines.substr(2));
        }
    }

    string metagenomic_reference_file = "../data/metagenomic_reference.fasta";
    if (filesystem::exists(metagenomic_reference_file)) {
        cout << "Metagenomic reference file already exists. Removing it..." << endl;
        filesystem::remove(metagenomic_reference_file);
    }
    create_metagenomic_reference(reference_files, metagenomic_reference_file);

    map<string,unordered_map<unsigned int, unsigned int>> ref_distributions;
    map<string,unordered_map<unsigned int, unsigned int>> frag_distributions;
    
    vector<analysis::Sequence> references = analysis::LoadSequences("../data/metagenomic_reference.fasta");

    for (const auto& ref : references) {
        vector<tuple<unsigned int, unsigned int, bool>> ref_minimizers = analysis::Minimize(ref.seq.c_str(), ref.seq.size(), k, w); 
        unordered_map<unsigned int, unsigned int> distribution_vector_ref_temp = distribution(ref_minimizers);

        cout << "Distribution for reference: " << ref.name << endl;

        for (const auto& [kmer, count] : distribution_vector_ref_temp) {
            cout << kmer << ": " << count << endl;
        }

        ref_distributions[ref.name] = distribution_vector_ref_temp;

    }

    for (const auto& frag : fragment_files) {
        vector<analysis::Sequence> fragments = analysis::LoadSequences("../data/fastq/" + frag);
        for (const auto& fragment : fragments) {
            vector<tuple<unsigned int, unsigned int, bool>> frag_minimizers = analysis::Minimize(fragment.seq.c_str(), fragment.seq.size(), k, w); 
            unordered_map<unsigned int, unsigned int> distribution_vector_frag_temp = distribution(frag_minimizers);

            cout << "Distribution for fragment: " << fragment.name << endl;

            for (const auto& [kmer, count] : distribution_vector_frag_temp) {
                cout << kmer << ": " << count << endl;
            }

        frag_distributions[fragment.name] = distribution_vector_frag_temp;

        }

    }

    for (auto& [frag_name, frag_dist] : frag_distributions) {
        double max_similarity = 0.0;
        string most_similar_reference;
        for (auto& [ref_name, ref_dist] : ref_distributions) {
            double similarity = cosine_similarity(ref_dist, frag_dist);
            cout << "Cosine similarity between " << ref_name << " and " << frag_name << ": " << similarity << endl;
            if (similarity > max_similarity) {
                max_similarity = similarity;
                most_similar_reference = ref_name;
            }
        }
        cout << "Most similar reference for fragment " << frag_name << ": " << most_similar_reference << " with similarity " << max_similarity << endl;
    }

    string reference_csv = "../data/reference_data.csv";
    string fragments_csv = "../data/fragment_data.csv";
    if (filesystem::exists(reference_csv)) {
        cout << "Reference CSV file already exists. Removing it..." << endl;
        filesystem::remove(reference_csv);
    }
    if (filesystem::exists(fragments_csv)) {
        cout << "Fragments CSV file already exists. Removing it..." << endl;
        filesystem::remove(fragments_csv);
    }

    export_to_csv(ref_distributions, reference_csv, k);
    export_to_csv(frag_distributions, fragments_csv, k);

    return 0;
}