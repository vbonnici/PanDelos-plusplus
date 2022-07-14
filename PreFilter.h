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
#include <omp.h>
#include <sys/time.h>
#include <cassert>
#include <typeinfo>
#include "include/global_options.h"

struct timeval tempo{};
double t1, t2, sum;


class PreFilter {
public:

    explicit PreFilter(const std::vector<std::string>& sequences,
                       const std::vector<std::vector<int>>& genome_sequencesid,
                       const int sequences_type,
                       std::ofstream* log_stream) :
        jaccard_threshold(0.8), sequences_type(sequences_type), kmer_size(6) {
        this->log_stream = log_stream;

        this->sequences = &sequences;
        this->genome_sequencesid = &genome_sequencesid;
        this->best_hits.reserve(sequences.size()*sequences.size());
        *this->log_stream << "byte riservati per prefilter best hits: " << sequences.size() * sequences.size() << std::endl;
    }

    void init_sequences_kmers() {
        this->sequences_kmers.reserve(this->sequences->size());
        for(int i = 0; i < this->sequences->size(); ++i) {
            std::array<unsigned int, 4095> array{0};

            this->sequences_kmers.emplace_back(array);
        }
    }

    void calculate_kmer_multiplicity() {
        if(this->sequences_type == 0)
            this->calculate_kmer_multiplicity_aminoacids();
        else
            this->calculate_kmer_multiplicity_nucleotides();

        *this->log_stream << "1 - kmer_multiplicity calcolate" << std::endl;
    }
    
    void calculate_best_hits() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;
        int index, i;

        if(!debug)
            omp_set_num_threads(omp_get_num_procs());

        for(index = 0; index < this->genome_sequencesid->size(); index++) {

            gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);
            std::vector<int> genome_a = this->genome_sequencesid->operator[](index);

            #pragma omp parallel for shared(genome_a) private(jaccard_similarity, counter_min, counter_max, sequence_a, sequence_b, value_a, value_b)
            for(i = index + 1; i < this->genome_sequencesid->size(); i++) {

                std::vector<std::pair<int, int>> best_hits_local;
                std::vector<int> genome_b = this->genome_sequencesid->operator[](i);

                best_hits_local.reserve(genome_a.size() * genome_b.size());

                for(auto &geneid_a: genome_a)
                    for(auto &geneid_b : genome_b) {
                        sequence_a = this->sequences->operator[](geneid_a);
                        sequence_b = this->sequences->operator[](geneid_b);

                        counter_min = 0;
                        counter_max = 0;
                        jaccard_similarity = 0;

                        std::array<unsigned int, 4095> kmer_array_a = this->sequences_kmers.operator[](geneid_a);
                        std::array<unsigned int, 4095> kmer_array_b = this->sequences_kmers.operator[](geneid_b);

                        if (!PreFilter::check_constraint(sequence_a, sequence_b))
                            continue;

                        for (int j = 0; j < 4095; ++j) {
                            value_a = kmer_array_a[j];
                            value_b = kmer_array_b[j];

                            if(value_a > value_b) {
                                counter_min += value_b;
                                counter_max += value_a;
                            }
                            else {
                                counter_min += value_a;
                                counter_max += value_b;
                            }
                        }

                        jaccard_similarity = (double) counter_min / counter_max;

                        if (counter_max > 0 && jaccard_similarity > this->jaccard_threshold) {
                            //*this->log_stream << geneid_a << " " << geneid_b << " prefilter jaccard similarity " << jaccard_similarity << std::endl;

                            best_hits_local.emplace_back(std::make_pair(geneid_a, geneid_b));
                        }
                    }

                #pragma omp critical
                this->best_hits.insert(this->best_hits.end(), std::make_move_iterator(best_hits_local.begin()), std::make_move_iterator(best_hits_local.end()));

            }

            gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);
            sum = (t2-t1);

            *this->log_stream << "Tempo computazione genoma " << index << " con gli altri genomi " << sum << std::endl;
            *this->log_stream << "Best hits in array: " << this->best_hits.size() << std::endl;
        }

        *this->log_stream << "1 - best hits calcolati " << std::endl;
    }

    std::vector<std::array<unsigned int, 4095>>& get_sequences_kmers() {
        return this->sequences_kmers;
    }

    std::vector<std::pair<int, int>>& get_best_hits() {
        return this->best_hits;
    }


private:
    const double jaccard_threshold;
    const int sequences_type; //0 amino acids, 1 nucleotides
    const int kmer_size;
    const std::vector<std::string>* sequences;
    const std::vector<std::vector<int>>* genome_sequencesid;
    std::vector<std::array<unsigned int, 4095>> sequences_kmers;
    std::vector<std::pair<int, int>> best_hits;
    std::ofstream* log_stream;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos;
    }

    [[nodiscard]] static bool aminoacid_is_valid(const std::string &str) {
        return str.find_first_not_of("FLIMVSPTAY*HQNKDECWRG") == std::string::npos;
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

                if(!aminoacid_is_valid(aminoacid))
                    continue;

                std::string kmer = PreFilter::aminoacid_to_nucleotides(aminoacid);

                int kmer_in_int = PreFilter::kmer_to_int(kmer);

                this->sequences_kmers.operator[](k).operator[](kmer_in_int) += 1;
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
