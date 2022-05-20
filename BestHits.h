#ifndef PANDELOS_PLUSPLUS_BESTHITS_H
#define PANDELOS_PLUSPLUS_BESTHITS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <array>

int kmer_array[4095] = {0};

class BestHits {
public:

    explicit BestHits(const std::vector<std::string>& sequences, const int flag) : jaccard_threshold(0.8), flag(flag) {
        this->sequences = &sequences;
    }

    void init_map_sequences_kmers() {
        for(auto &sequence: *this->sequences) {
            this->map_sequences_kmers.insert(std::make_pair(sequence, kmer_array));
        }
    }

    void calculate_kmer_multiplicity_nucleotides() {
        int counter = 0;
        for(auto &i : this->map_sequences_kmers) {
            std::vector<std::string> sequence = split_string(i.first, '\n');

            for(int window = 0; window < sequence[1].length(); ++window) {
                std::string kmer = sequence[1].substr(window, 6);

                if(kmer_is_valid(kmer)) {
                    //std::cout << kmer << std::endl;
                    int kmer_in_int = kmer_to_int(kmer);
                    //std::cout << kmer_in_int << std::endl;

                    auto* index = i.second + kmer_in_int;
                    *index += 1;

                    std::cout << *index << std::endl;
                }
            }
            ++counter;
            std::cout << sequence[0] << counter << std::endl;
        }
    }

    std::unordered_map<std::string, int*>& get_map_sequences_kmers() {
        return this->map_sequences_kmers;
    }

private:
    const double jaccard_threshold;
    const std::vector<std::string>* sequences;
    const int flag;
    std::unordered_map<std::string, int*> map_sequences_kmers;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> map_best_hits;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == 6 && str.find_first_not_of("ACGT") == std::string::npos;
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
};
#endif //PANDELOS_PLUSPLUS_BESTHITS_H
