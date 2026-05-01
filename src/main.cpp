#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <fstream>
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


int main(int argc, char *argv[]){
    cout << "Starting program..." << endl;

    string configuration = "config.txt";
    string lines;
    unsigned int k = 10;
    unsigned int w = 3;
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
    
    vector<analysis::Sequence> references = analysis::LoadSequences("../data/metagenomic_reference.fasta");
    
    for (const auto& frag : fragment_files) {
        vector<analysis::Sequence> fragments = analysis::LoadSequences("../data/fastq/" + frag);
    }

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