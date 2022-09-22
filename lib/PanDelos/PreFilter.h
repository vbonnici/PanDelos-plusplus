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
#include "../../conf/conf.h"
#include "../../include/Helper.h"

struct timeval chrono{};
double t1, t2, sum;


class PreFilter {
public:

    explicit PreFilter(const std::vector<std::string>& sequences,
                       const std::vector<std::vector<int>>& genome_sequencesid,
                       const int sequences_type,
                       std::ofstream* log_stream) :
        jaccard_threshold(0.7), sequences_type(sequences_type), kmer_size(6) {
        this->log_stream = log_stream;

        this->sequences = &sequences;
        this->genome_sequencesid = &genome_sequencesid;
    }

    void init_sequences_kmers() {
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
    }

    void find_candidate_sequences(const std::vector<std::pair<int, int>>& gene_id_pair_in, double jaccard_threshold_in) {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;

        omp_set_num_threads(omp_get_num_procs());

        #pragma omp parallel
        {
            std::vector<std::pair<int, int>> candidate_sequences_local;
            #pragma omp for private(jaccard_similarity, counter_min, counter_max, sequence_a, sequence_b, value_a, value_b) schedule(static)
            for(auto &i : gene_id_pair_in) {

                sequence_a = this->sequences->operator[](i.first);
                sequence_b = this->sequences->operator[](i.second);

                counter_min = 0;
                counter_max = 0;
                jaccard_similarity = 0;

                std::array<unsigned int, 4095> kmer_array_a = this->sequences_kmers.operator[](i.first);
                std::array<unsigned int, 4095> kmer_array_b = this->sequences_kmers.operator[](i.second);

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

                if (counter_max > 0 && jaccard_similarity > jaccard_threshold_in)
                    candidate_sequences_local.emplace_back(std::make_pair(i.first, i.second));
            }

            for (int t = 0; t < omp_get_num_threads(); t++) {
                #pragma omp barrier
                if (t == omp_get_thread_num()) {
                    this->candidate_sequences.insert(this->candidate_sequences.end(), std::make_move_iterator(candidate_sequences_local.begin()), std::make_move_iterator(candidate_sequences_local.end()));
                }
            }
        };

    }
    
    void find_candidate_sequences() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;
        int index, i;

        omp_set_num_threads(omp_get_num_procs());

        for(index = 0; index < this->genome_sequencesid->size(); index++) {

            gettimeofday(&chrono,nullptr); t1 = chrono.tv_sec+(chrono.tv_usec/1000000.0);
            std::vector<int> genome_a = this->genome_sequencesid->operator[](index);

            #pragma omp parallel for shared(genome_a) private(jaccard_similarity, counter_min, counter_max, sequence_a, sequence_b, value_a, value_b)
            for(i = index + 1; i < this->genome_sequencesid->size(); i++) {

                std::vector<std::pair<int, int>> candidate_sequences_local;
                std::vector<int> genome_b = this->genome_sequencesid->operator[](i);

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

                        if (counter_max > 0 && jaccard_similarity > this->jaccard_threshold)
                            candidate_sequences_local.emplace_back(std::make_pair(geneid_a, geneid_b));
                    }

                #pragma omp critical
                this->candidate_sequences.insert(this->candidate_sequences.end(), std::make_move_iterator(candidate_sequences_local.begin()), std::make_move_iterator(candidate_sequences_local.end()));

            }

            gettimeofday(&chrono,nullptr); t2 = chrono.tv_sec+(chrono.tv_usec/1000000.0);
            sum = (t2-t1);

            //*this->log_stream << "Genome computation time " << index << " with other genomes " << sum << std::endl;
            //*this->log_stream << "Candidate sequences in array: " << this->candidate_sequences.size() << std::endl;
        }
    }

    std::vector<std::array<unsigned int, 4095>>& get_sequences_kmers() {
        return this->sequences_kmers;
    }

    std::vector<std::pair<int, int>>& get_candidate_sequences() {
        return this->candidate_sequences;
    }

private:
    const double jaccard_threshold;
    const int sequences_type;
    const int kmer_size;
    const std::vector<std::string>* sequences;
    const std::vector<std::vector<int>>* genome_sequencesid;
    std::vector<std::array<unsigned int, 4095>> sequences_kmers;
    std::vector<std::pair<int, int>> candidate_sequences;
    std::ofstream* log_stream;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos;
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
        std::string sequence;
        std::string kmer;
        int kmer_in_int;

        #pragma omp parallel for private(sequence, kmer, kmer_in_int)
        for(int k = 0; k < this->sequences_kmers.size(); ++k) {
            sequence = this->sequences->operator[](k);

            for(int window = 0; window < sequence.length() - this->kmer_size + 1; window++) {
                kmer = sequence.substr(window, this->kmer_size);

                if(this->kmer_is_valid(kmer)) {
                    kmer_in_int = PreFilter::kmer_to_int(kmer);

                    this->sequences_kmers.operator[](k).operator[](kmer_in_int) += 1;
                }
            }
        }
    }

    void calculate_kmer_multiplicity_aminoacids() {
        std::string sequence;
        std::string aminoacid;
        std::string kmer;
        int kmer_in_int;

        #pragma omp parallel for private(sequence, aminoacid, kmer, kmer_in_int)
        for(int k = 0; k < this->sequences_kmers.size(); ++k) {
            sequence = this->sequences->operator[](k);

            for(int window = 0; window < sequence.length() - this->kmer_size + 1; window++) {
                aminoacid = sequence.substr(window, 2);

                kmer = Helper::aminoacid_to_nucleotides(aminoacid);

                kmer_in_int = PreFilter::kmer_to_int(kmer);

                this->sequences_kmers.operator[](k).operator[](kmer_in_int) += 1;
            }
        }
    }
};
#endif //PANDELOS_PLUSPLUS_PREFILTER_H
