//
// Created by Giandonato Inverso on 06/07/2022.
//

#ifndef PANDELOS_PLUSPLUS_PARALOG_H
#define PANDELOS_PLUSPLUS_PARALOG_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <algorithm>
#include "include/kvalue.h"
#include <omp.h>
#include <sys/time.h>
#include <array>

class Paralog {
public:
    explicit Paralog(const std::vector<std::string>& sequences_prefilter,
                     std::vector<std::vector<int>>& genome_sequencesid,
                     const int sequences_type,
                     const int kmer_size,
                     std::vector<std::tuple<int, int, double>>& vector_tuple_bbh,
                     std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size){

        this->log_stream = log_stream;
        this->sequences = &sequences_prefilter;
        this->genome_sequencesid = &genome_sequencesid;
        this->vector_tuple_bbh = &vector_tuple_bbh;
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

        *this->log_stream << "3 - kmer_multiplicity calcolate" << std::endl;
    }

    void calculate_paralog() {
        this->genome_minimum_jaccard.reserve(this->genome_sequencesid->size());

        for(int i = 0; i < this->genome_sequencesid->size(); ++i)
            this->genome_minimum_jaccard.push_back(1.0);

        this->compute_minimum_jaccard();
        this->calculate_paralog_best_hits();
    }

    std::vector<std::tuple<int, int, double>>& get_paralog_best_hits() {
        return this->paralog_best_hits;
    }

private:
    std::ofstream* log_stream;
    int kmer_size;
    const int sequences_type; //0 amino acids, 1 nucleotides
    std::vector<std::vector<int>>* genome_sequencesid;
    std::vector<double> genome_minimum_jaccard;
    std::vector<std::array<unsigned int, 4095>> sequences_kmers;
    const std::vector<std::string>* sequences;
    std::vector<std::tuple<int, int, double>>* vector_tuple_bbh;
    std::vector<std::tuple<int, int, double>> paralog_best_hits;


    void compute_minimum_jaccard() {
        std::vector<int> limit_genes_id;

        for(int i = 0; i < this->genome_sequencesid->size(); ++i) {
            auto genes = this->genome_sequencesid->operator[](i);

            auto max = max_element(std::begin(genes), std::end(genes));

            limit_genes_id.push_back(*max);
        }

        for(int i = 0; i < this->genome_minimum_jaccard.size(); ++i) {
            for(auto &tuple : *this->vector_tuple_bbh) {
                if(std::get<0>(tuple) <= limit_genes_id.operator[](i))
                    if(std::get<2>(tuple) < this->genome_minimum_jaccard.operator[](i))
                        this->genome_minimum_jaccard.operator[](i) = std::get<2>(tuple);
            }
        }
    }

    void calculate_paralog_best_hits() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;
        int index, i;

        omp_set_num_threads(omp_get_num_procs());

        #pragma omp parallel for private(index, i, t1, t2, sum, jaccard_similarity, counter_min, counter_max, sequence_a, sequence_b, value_a, value_b)
        for(index = 0; index < this->genome_sequencesid->size(); index++) {

            gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);
            std::vector<int> genome_a = this->genome_sequencesid->operator[](index);

            //#pragma omp parallel for shared(genome_a) private(jaccard_similarity, counter_min, counter_max, sequence_a, sequence_b, value_a, value_b)

            std::vector<std::tuple<int, int, double>> best_hits_local;
            std::vector<int> genome_b = this->genome_sequencesid->operator[](index);

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

                    if (!Paralog::check_constraint(sequence_a, sequence_b))
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

                    if (counter_max > 0 && jaccard_similarity > this->genome_minimum_jaccard.operator[](index)) {
                        *this->log_stream << geneid_a << " " << geneid_b << " paralog jaccard similarity " << jaccard_similarity << std::endl;

                        best_hits_local.emplace_back(std::make_tuple(geneid_a, geneid_b, jaccard_similarity));
                    }
                }

            #pragma omp critical
            this->paralog_best_hits.insert(this->paralog_best_hits.end(), std::make_move_iterator(best_hits_local.begin()), std::make_move_iterator(best_hits_local.end()));

            gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);
            sum = (t2-t1);

            *this->log_stream << "Tempo computazione genoma " << index << " con se stesso " << sum << std::endl;
            *this->log_stream << "Paralog Best hits in array: " << this->paralog_best_hits.size() << std::endl;
        }

        *this->log_stream << "3 - paralog best hits calcolati " << std::endl;
    }

    /*TODO aggiungere all'helper */

    static bool check_constraint(std::string& sequence_a, std::string& sequence_b) {

        if(sequence_a.length() >= sequence_b.length()*2 || sequence_b.length() >= sequence_a.length()*2)
            return false;

        return true;
    }

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

    void calculate_kmer_multiplicity_nucleotides() {
        for(int k = 0; k < this->sequences_kmers.size(); ++k) {

            std::string sequence = this->sequences->operator[](k);

            for(int window = 0; window < sequence.length() - this->kmer_size + 1; window++) {
                std::string kmer = sequence.substr(window, this->kmer_size);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = Paralog::kmer_to_int(kmer);

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
                std::string kmer = Paralog::aminoacid_to_nucleotides(aminoacid);

                if(kmer_is_valid(kmer)) {
                    int kmer_in_int = Paralog::kmer_to_int(kmer);

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


#endif //PANDELOS_PLUSPLUS_PARALOG_H
