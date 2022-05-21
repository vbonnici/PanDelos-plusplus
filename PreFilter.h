#ifndef PANDELOS_PLUSPLUS_PREFILTER_H
#define PANDELOS_PLUSPLUS_PREFILTER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <array>

class PreFilter {
public:

    explicit PreFilter(const std::vector<std::string>& sequences, const int flag) :
        jaccard_threshold(0.8), flag(flag), kmer_size(6) {
        this->sequences = &sequences;
    }

    void init_map_sequences_kmers() {
        for(auto &sequence: *this->sequences) {
            std::array<unsigned int, 4095> array{0};

            this->map_sequences_kmers.insert(std::make_pair(sequence, array));
        }
    }

    void calculate_kmer_multiplicity() {
        if(this->flag == 0)
            this->calculate_kmer_multiplicity_aminoacids();
        else
            this->calculate_kmer_multiplicity_nucleotides();
    }

    void calculate_best_hits() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;

        for(auto &a : this->map_sequences_kmers) {
            for(auto &b : this->map_sequences_kmers) {
                std::vector<std::string> sequence_a = PreFilter::split_string(a.first, '\n');
                std::vector<std::string> sequence_b = PreFilter::split_string(b.first, '\n');

                counter_min = 0;
                counter_max = 0;

                if (!PreFilter::check_constraint(sequence_a, sequence_b))
                    continue;

                for (int j = 0; j < 4095; ++j) {

                    value_a = a.second[j];
                    value_b = b.second[j];

                    //std::cout << "sequence_a " << sequence_a[0] << " sequence_b " << sequence_b[0] << std::endl;
                    //std::cout << "value_a " << value_a << " value_b " << value_b << std::endl;

                    counter_min += std::min(value_a, value_b);
                    counter_max += std::max(value_a, value_b);

                    //std::cout << "min " << counter_min << " max " << counter_max << std::endl;

                }

                jaccard_similarity = (double) counter_min / counter_max;
                std::cout << "jaccard similarity " << jaccard_similarity << std::endl;

                if (isfinite(jaccard_similarity) && jaccard_similarity > this->jaccard_threshold) {
                    std::unordered_map<std::string, double> temp;
                    temp.insert(std::make_pair(sequence_b[0], jaccard_similarity));

                    this->map_best_hits.insert(std::make_pair(sequence_a[0], temp));
                }
            }
        }
        std::cout << "best hits calcolati " << std::endl;
    }

    std::unordered_map<std::string, std::array<unsigned int, 4095>>& get_map_sequences_kmers() {
        return this->map_sequences_kmers;
    }

    std::unordered_map<std::string, std::unordered_map<std::string, double>>& get_map_best_hits() {
        return this->map_best_hits;
    }

private:
    const double jaccard_threshold;
    const int flag; //0 amino acids, 1 nucleotides
    const int kmer_size;
    const std::vector<std::string>* sequences;
    std::unordered_map<std::string, std::array<unsigned int, 4095>> map_sequences_kmers;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> map_best_hits;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos;
    }

    static std::vector<std::string> split_string(std::string const &str, const char delim) {
        std::vector<std::string> out;
        size_t start;
        size_t end = 0;

        while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
            end = str.find(delim, start);
            out.push_back(str.substr(start, end - start));
        }

        return out;
    }

    static int kmer_to_int(std::string& kmer) {
        int bitmap = 0;
        int A = 0b00;
        int C = 0b01;
        int G = 0b10;
        int T = 0b11;

        for(int a = 0; a < kmer.length(); ++a) {

            if(reinterpret_cast<char>(kmer[a]) == 'A')
                bitmap = (bitmap << 2) | A;

            else if(reinterpret_cast<char>(kmer[a]) == 'C')
                bitmap = (bitmap << 2) | C;

            else if(reinterpret_cast<char>(kmer[a]) == 'G')
                bitmap = (bitmap << 2) | G;

            else if (reinterpret_cast<char>(kmer[a]) == 'T')
                bitmap = (bitmap << 2) | T;
        }

        return bitmap;
    }

    static bool check_constraint(std::vector<std::string>& sequence_a, std::vector<std::string>& sequence_b) {
        if(sequence_a[0] == sequence_b[0])
            return false;

        if(sequence_a[0].length() >= sequence_b[0].length()*2 || sequence_b[0].length() >= sequence_a[0].length()*2)
            return false;

        return true;
    }

    void calculate_kmer_multiplicity_nucleotides() {
        for(auto &i : this->map_sequences_kmers) {
            std::vector<std::string> sequence = PreFilter::split_string(i.first, '\n');

            for(int window = 0; window < sequence[1].length(); ++window) {
                std::string kmer = sequence[1].substr(window, this->kmer_size);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = PreFilter::kmer_to_int(kmer);

                    i.second[kmer_in_int] += 1;
                }
            }
        }
    }

    void calculate_kmer_multiplicity_aminoacids() {
        for(auto &i : this->map_sequences_kmers) {
            std::vector<std::string> sequence = PreFilter::split_string(i.first, '\n');

            for(int window = 0; window < sequence[1].length(); ++window) {
                std::string aminoacid = sequence[1].substr(window, 2);
                std::string kmer = PreFilter::aminoacid_to_nucleotides(aminoacid);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = PreFilter::kmer_to_int(kmer);

                    i.second[kmer_in_int] += 1;
                }
            }
        }
    }

    static std::string aminoacid_to_nucleotides(std::string& aminoacid) {
        std::string kmer;

        for(int a = 0; a < aminoacid.length(); ++a) {

            if(reinterpret_cast<char>(aminoacid[a]) == 'F')
                kmer += "TTT";

            else if(reinterpret_cast<char>(aminoacid[a]) == 'L')
                kmer += "TTA";

            else if(reinterpret_cast<char>(aminoacid[a]) == 'I')
                kmer += "ATT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'M')
                kmer += "ATG";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'V')
                kmer += "GTT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'S')
                kmer += "TCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'P')
                kmer += "CCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'T')
                kmer += "ACT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'A')
                kmer += "GCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'Y')
                kmer += "TAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == '*')
                kmer += "TAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'H')
                kmer += "CAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'Q')
                kmer += "CAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'N')
                kmer += "AAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'K')
                kmer += "AAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'D')
                kmer += "GAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'E')
                kmer += "GAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'C')
                kmer += "TGT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'W')
                kmer += "TGG";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'R')
                kmer += "CGT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'G')
                kmer += "GGG";
        }

        return kmer;
    }
};
#endif //PANDELOS_PLUSPLUS_PREFILTER_H
