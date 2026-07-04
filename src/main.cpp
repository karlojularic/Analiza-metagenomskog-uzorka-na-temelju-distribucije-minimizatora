#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <map>
#include "loading_sequences.hpp"
#include "minimizers.hpp"

using namespace std;

auto distribution(const vector<pair<unsigned int, unsigned int>>& minimizers) {
    unordered_map<unsigned int, unsigned int> distribution_vector;

    for (const auto& [kmer, pos] : minimizers) {
        distribution_vector[kmer]++;
    }

    return distribution_vector;

}

void create_metagenomic_reference(const vector<string>& reference_files, const string& output_file) {
    ofstream out(output_file);

    for (const auto& file : reference_files) {
        ifstream in("../data/reference_files/" + file);
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
            dot += (double)count * (double)frag.at(kmer);
        }
        denom_a += (double)count * (double)count;
    }

    for (const auto& [kmer, count] : frag) {
        denom_b += (double)count * (double)count;
    }

    return dot / (sqrt(denom_a) * sqrt(denom_b));

}

void export_to_csv(const map<string, unordered_map<unsigned int, unsigned int>>& distributions, const string& output_file, unsigned int k) {
    ofstream out(output_file);

    for (const auto& [name, dist] : distributions) {
        out << name;
        if (&name != &distributions.rbegin()->first) {
            out << ",";
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
                out << ",";
            }
        }
        out << endl;
    }
    out.close();
}

void export_classification_to_csv(const unordered_map<string, pair<string, double>>& classification_map, const unordered_map<string, string>& fragment_to_reference_map, const string& output_file) {
    ofstream out(output_file);
    out << "Fragment,Most Similar Reference,Actual Reference,Cosine Similarity" << endl;
    for (const auto& [fragment, references] : classification_map) {
        const auto& [reference, similarity] = references;
        out << fragment << "," << reference << "," << fragment_to_reference_map.at(fragment) << "," << similarity << endl;
    }
    out.close();
}


int main(int argc, char *argv[]){
    cout << "Starting program..." << endl;

    string configuration = "../config.txt";
    string line;
    unsigned int k = 6;
    unsigned int w = 10;
    vector<string> reference_files;
    vector<string> fragment_files;

    if(!(filesystem::exists(configuration))) {
        cerr << "Configuration file does not exist: " << configuration << endl;
        return 1;
    }

    ifstream config_file(configuration);
    while (getline(config_file, line)) {
        if (line.rfind("reference_files=", 0) == 0) {
            string files_str = line.substr(16);
            size_t pos = 0;
            while ((pos = files_str.find(',')) != string::npos) {
                reference_files.push_back(files_str.substr(0, pos));
                files_str.erase(0, pos + 1);
            }
            reference_files.push_back(files_str);
        } else if (line.rfind("fragment_files=", 0) == 0) {
            string files_str = line.substr(15);
            size_t pos = 0;
            while ((pos = files_str.find(',')) != string::npos) {
                fragment_files.push_back(files_str.substr(0, pos));
                files_str.erase(0, pos + 1);
            }
            fragment_files.push_back(files_str);
        } else if (line.rfind("k=", 0) == 0) {
            k = stoi(line.substr(2));
        } else if (line.rfind("w=", 0) == 0) {
            w = stoi(line.substr(2));
        }
    }

    string metagenomic_reference_file = "../output/metagenomic_reference.fasta";

    if (!filesystem::exists("../output")) {
        cout << "Output directory does not exist. Creating it..." << endl;
        filesystem::create_directory("../output");
    }
    if (filesystem::exists(metagenomic_reference_file)) {
        cout << "Metagenomic reference file already exists. Removing it..." << endl;
        filesystem::remove(metagenomic_reference_file);
    }
    create_metagenomic_reference(reference_files, metagenomic_reference_file);

    map<string,unordered_map<unsigned int, unsigned int>> ref_fwd_distributions;
    map<string,unordered_map<unsigned int, unsigned int>> ref_rc_distributions;
    map<string,unordered_map<unsigned int, unsigned int>> frag_distributions;
    
    vector<analysis::Sequence> references = analysis::LoadSequences("../output/metagenomic_reference.fasta");

    for (const auto& ref : references) {
        vector<pair<unsigned int, unsigned int>> ref_fwd_minimizers = analysis::Minimize(ref.seq.c_str(), ref.seq.size(), k, w); 
        unordered_map<unsigned int, unsigned int> ref_fwd_dist_temp = distribution(ref_fwd_minimizers);

        cout << "Distribution for reference: " << ref.name << endl;

        for (const auto& [kmer, count] : ref_fwd_dist_temp) {
            cout << kmer << ": " << count << endl;
        }

        string ref_name = ref.name.substr(0, ref.name.find(' '));
        ref_fwd_distributions[ref_name] = ref_fwd_dist_temp;

        string ref_rc_seq = analysis::getReverseChain(ref.seq);
        vector<pair<unsigned int, unsigned int>> ref_rc_minimizers = analysis::Minimize(ref_rc_seq.c_str(), ref_rc_seq.size(), k, w);
        unordered_map<unsigned int, unsigned int> ref_rc_dist_temp = distribution(ref_rc_minimizers);
        ref_rc_distributions[ref_name] = ref_rc_dist_temp;

    }

    unordered_map<string, string> fragment_to_reference_map;

    for (const auto& frag : fragment_files) {
        vector<analysis::Sequence> fragments = analysis::LoadSequences("../data/fragments_files/" + frag);
        for (const auto& fragment : fragments) {
            vector<pair<unsigned int, unsigned int>> frag_minimizers = analysis::Minimize(fragment.seq.c_str(), fragment.seq.size(), k, w); 
            unordered_map<unsigned int, unsigned int> frag_dist_temp = distribution(frag_minimizers);

            cout << "Distribution for fragment: " << fragment.name << endl;

            for (const auto& [kmer, count] : frag_dist_temp) {
                cout << kmer << ": " << count << endl;
            }

            frag_distributions[fragment.name] = frag_dist_temp;
            fragment_to_reference_map[fragment.name] = frag.substr(0, frag.find("_trimmed.fastq"));

        }

    }

    unordered_map<string, pair<string, double>> best_match_map_fwd;
    for (auto& [frag_name, frag_dist] : frag_distributions) {
        double max_similarity = 0.0;
        string most_similar_reference;
        for (auto& [ref_name, ref_dist] : ref_fwd_distributions) {
            double similarity = cosine_similarity(ref_dist, frag_dist);
            cout << "Cosine similarity between " << ref_name << " and " << frag_name << ": " << similarity << endl;
            if (similarity > max_similarity) {
                max_similarity = similarity;
                most_similar_reference = ref_name;
            }
        }
        best_match_map_fwd[frag_name] = {most_similar_reference, max_similarity};
        cout << "Most similar reference for fragment " << frag_name << ": " << most_similar_reference << " with similarity " << max_similarity << endl;
    }

    unordered_map<string, pair<string, double>> best_match_map_rc;
    for (auto& [frag_name, frag_dist] : frag_distributions) {
        double max_similarity = 0.0;
        string most_similar_reference;
        for (auto& [ref_name, ref_dist] : ref_rc_distributions) {
            double similarity = cosine_similarity(ref_dist, frag_dist);
            cout << "Cosine similarity between " << ref_name << " and " << frag_name << ": " << similarity << endl;
            if (similarity > max_similarity) {
                max_similarity = similarity;
                most_similar_reference = ref_name;
            }
        }
        best_match_map_rc[frag_name] = {most_similar_reference, max_similarity};
        cout << "Most similar reference for fragment " << frag_name << ": " << most_similar_reference << " with similarity " << max_similarity << endl;
    }

    unordered_map<string, pair<string, double>> combined_best_match_map;
    for (const auto& [frag_name_fwd, frag_match_fwd] : best_match_map_fwd) {
        auto matching_rc = best_match_map_rc.at(frag_name_fwd);
        if (frag_match_fwd.second >= matching_rc.second) {
            combined_best_match_map[frag_name_fwd] = frag_match_fwd;
        } else {
            combined_best_match_map[frag_name_fwd] = matching_rc;
        }
    }

    string reference_fwd_csv = "../output/reference_data_fwd.csv";
    string reference_rc_csv = "../output/reference_data_rc.csv";
    string fragments_csv = "../output/fragment_data.csv";
    string classification_csv = "../output/classification_data.csv";
    if (filesystem::exists(reference_fwd_csv)) {
        cout << "Reference CSV file already exists. Removing it..." << endl;
        filesystem::remove(reference_fwd_csv);
    }
    if (filesystem::exists(reference_rc_csv)) {
        cout << "Reference RC CSV file already exists. Removing it..." << endl;
        filesystem::remove(reference_rc_csv);
    }
    if (filesystem::exists(fragments_csv)) {
        cout << "Fragments CSV file already exists. Removing it..." << endl;
        filesystem::remove(fragments_csv);
    }
    if(filesystem::exists(classification_csv)) {
        cout << "Classification CSV file already exists. Removing it..." << endl;
        filesystem::remove(classification_csv);
    }

    export_to_csv(ref_fwd_distributions, reference_fwd_csv, k);
    export_to_csv(ref_rc_distributions, reference_rc_csv, k);
    export_to_csv(frag_distributions, fragments_csv, k);
    export_classification_to_csv(combined_best_match_map, fragment_to_reference_map, classification_csv);
    cout << "\nRunning Python analysis script..." << endl;
    system("cd .. && venv/bin/python3 analysis.py");

    cout << "\nGenerating PDF report..." << endl;
    system("cd .. && Rscript -e \"rmarkdown::render('distribution_analysis.Rmd', output_file='output/distribution_report.pdf')\"");

    return 0;
}