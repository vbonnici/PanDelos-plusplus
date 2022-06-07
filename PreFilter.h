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
#include <algorithm>

class PreFilter {
public:

    explicit PreFilter(const std::vector<std::string>& sequences, const std::vector<std::vector<int>>& genome_sequencesid, const int flag) :
        jaccard_threshold(0.8), flag(flag), kmer_size(6) {
        this->sequences = &sequences;
        this->genome_sequencesid = &genome_sequencesid;
    }

    void init_sequences_kmers() {
        this->sequences_kmers.reserve(this->sequences->size());
        for(int i = 0; i < this->sequences->size(); ++i) {
            std::array<unsigned int, 4095> array{0};

            this->sequences_kmers.emplace_back(array);
        }
    }

    void calculate_kmer_multiplicity() {
        if(this->flag == 0)
            this->calculate_kmer_multiplicity_aminoacids();
        else
            this->calculate_kmer_multiplicity_nucleotides();
    }

    void calculate_best_hits_v2() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;

        int index ;

        for(index = 0; index < this->genome_sequencesid->size(); ++index) {
            for(int i = index + 1; i < this->genome_sequencesid->size(); ++i) {

                std::vector<int> genome_a = this->genome_sequencesid->operator[](index);
                std::vector<int> genome_b = this->genome_sequencesid->operator[](i);

                for(auto & a: genome_a)
                    for(auto & b : genome_b) {

                        sequence_a = this->sequences->operator[](a);
                        sequence_b = this->sequences->operator[](b);

                        counter_min = 0;
                        counter_max = 0;

                        std::array<unsigned int, 4095> kmer_array_a = this->sequences_kmers.operator[](a);
                        std::array<unsigned int, 4095> kmer_array_b = this->sequences_kmers.operator[](b);

                        if (!PreFilter::check_constraint(sequence_a, sequence_b))
                            continue;

                        for (int j = 0; j < 4095; ++j) {

                            value_a = kmer_array_a[j];
                            value_b = kmer_array_b[j];

                            //std::cout << "value_a " << value_a << " value_b " << value_b << std::endl;

                            if(value_a > value_b) {
                                counter_min += value_b;
                                counter_max += value_a;
                            }
                            else {
                                counter_min += value_a;
                                counter_max += value_b;
                            }

                            //std::cout << "min " << counter_min << " max " << counter_max << std::endl;

                        }

                        jaccard_similarity = (double) counter_min / counter_max;
                        //std::cout << "prefilter jaccard similarity " << jaccard_similarity << std::endl;

                        if (counter_max > 0 && jaccard_similarity > this->jaccard_threshold) {
                            this->best_hits.emplace_back(std::make_pair(a, b));
                            //std::cout << a << " best hit con " << b << std::endl;
                        }
                    }
            }
        }

        std::cout << "prefilter best hits calcolati " << std::endl;
    }

    std::vector<std::array<unsigned int, 4095>>& get_sequences_kmers() {
        return this->sequences_kmers;
    }

    std::vector<std::pair<int, int>>& get_best_hits() {
        return this->best_hits;
    }

private:
    const double jaccard_threshold;
    const int flag; //0 amino acids, 1 nucleotides
    const int kmer_size;
    const std::vector<std::string>* sequences;
    const std::vector<std::vector<int>>* genome_sequencesid;
    std::vector<std::array<unsigned int, 4095>> sequences_kmers;
    std::vector<std::pair<int, int>> best_hits;

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

    static bool check_constraint(std::string& sequence_a, std::string& sequence_b) {
        if(sequence_a == sequence_b)
            return false;

        if(sequence_a.length() >= sequence_b.length()*2 || sequence_b.length() >= sequence_a.length()*2)
            return false;

        return true;
    }

    void calculate_kmer_multiplicity_nucleotides() {
        for(int k = 0; k < this->sequences_kmers.size(); ++k) {

            std::string sequence = this->sequences->operator[](k);

            for(int window = 0; window < sequence.length() - this->kmer_size + 1; window++) {
                std::string kmer = sequence.substr(window, this->kmer_size);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = PreFilter::kmer_to_int(kmer);

                    this->sequences_kmers.operator[](k).operator[](kmer_in_int) += 1;
                }
            }
        }
    }

    void calculate_kmer_multiplicity_aminoacids() {
        for(int k = 0; k < this->sequences_kmers.size(); ++k) {

            std::string sequence = this->sequences->operator[](k);

            for(int window = 0; window < sequence.length() - this->kmer_size + 1; window++) {
                std::string aminoacid = sequence.substr(window, 2);
                std::string kmer = PreFilter::aminoacid_to_nucleotides(aminoacid);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = PreFilter::kmer_to_int(kmer);

                    this->sequences_kmers.operator[](k).operator[](kmer_in_int) += 1;
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
