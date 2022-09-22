#ifndef PANDELOS_PLUSPLUS_PARALOGUES_H
#define PANDELOS_PLUSPLUS_PARALOGUES_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include "../../conf/conf.h"
#include "Homologues.h"
#include "PreFilter.h"

class Paralogues {
public:
    explicit Paralogues(std::vector<std::string> input_sequences,
                     std::vector<std::vector<int>>& genome_sequencesid,
                     std::vector<std::pair<int, int>>& genes_id_interval,
                     const int sequences_type,
                     int kmer_size,
                     std::vector<std::tuple<int, int, double>>& bidirectional_best_hits,
                     std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size) {

        this->log_stream = log_stream;
        this->sequences = std::move(input_sequences);
        this->genome_sequencesid = genome_sequencesid;
        this->bidirectional_best_hits = &bidirectional_best_hits;
        this->genes_id_interval = &genes_id_interval;
        this->genome_counter = this->genome_sequencesid.size();
    }


    void find_paralogues() {
        for(int i = 0; i < this->genome_counter; ++i)
            this->genome_minimum_jaccard.push_back(1.0);

        this->calculate_minimum_jaccard();
        this->find_paralogues_best_hits();
    }

    std::vector<std::tuple<int, int, double>>& get_paralogues_best_hits() {
        return this->paralogues_best_hits;
    }

private:
    std::ofstream* log_stream;
    int kmer_size;
    const int sequences_type;
    int genome_counter;
    std::vector<std::vector<int>> genome_sequencesid;
    std::vector<double> genome_minimum_jaccard;
    std::vector<std::string> sequences;
    std::vector<std::tuple<int, int, double>>* bidirectional_best_hits;
    std::vector<std::tuple<int, int, double>> paralogues_best_hits;
    std::vector<std::pair<int, int>>* genes_id_interval;


    void calculate_minimum_jaccard() {
        for(int i = 0; i < this->genome_counter; ++i) {
            for(auto &tuple : *this->bidirectional_best_hits) {
                if(std::get<0>(tuple) >= std::get<0>(this->genes_id_interval->operator[](i)) && std::get<0>(tuple) <= std::get<1>(this->genes_id_interval->operator[](i))) {
                    if(std::get<2>(tuple) < this->genome_minimum_jaccard.operator[](i))
                        this->genome_minimum_jaccard.operator[](i) = std::get<2>(tuple);
                }

                if(std::get<1>(tuple) >= std::get<0>(this->genes_id_interval->operator[](i)) && std::get<1>(tuple) <= std::get<1>(this->genes_id_interval->operator[](i))) {
                    if(std::get<2>(tuple) < this->genome_minimum_jaccard.operator[](i))
                        this->genome_minimum_jaccard.operator[](i) = std::get<2>(tuple);
                }
            }
        }

        /*for(int i = 0; i < this->genome_counter; ++i) {
            *this->log_stream << "genome id " << i << " minimum jaccard " << this->genome_minimum_jaccard.operator[](i) << std::endl;
        }*/
    }

    void find_paralogues_best_hits() {
        int min_gene_id;
        int max_gene_id;

        for(int i = 0; i < this->genome_counter; ++i) {
            min_gene_id = std::get<0>(this->genes_id_interval->operator[](i));
            max_gene_id = std::get<1>(this->genes_id_interval->operator[](i));

            std::vector<std::pair<int, int>> gene_id_pair;

            for(int a = min_gene_id; a <= max_gene_id; ++a)
                for(int b = a+1; b <= max_gene_id; ++b)
                    if(this->check_constraint(a, b))
                        gene_id_pair.emplace_back(std::make_pair(a, b));


            PreFilter pre_filter = PreFilter(this->sequences, this->genome_sequencesid, this->sequences_type, this->log_stream);
            pre_filter.init_sequences_kmers();
            pre_filter.calculate_kmer_multiplicity();
            pre_filter.find_candidate_sequences(gene_id_pair, this->genome_minimum_jaccard.operator[](i));
            auto pre_filter_candidate_sequences = pre_filter.get_candidate_sequences();

            Homologues homologues = Homologues(this->sequences, pre_filter_candidate_sequences, this->sequences_type, this->kmer_size, this->log_stream);

            homologues.init_sequences_kmers();
            homologues.calculate_kmer_multiplicity();
            homologues.find_candidate_sequences(this->genome_minimum_jaccard.operator[](i));

            auto homologues_candidate_sequences = homologues.get_candidate_sequences();

            for(auto &it : homologues_candidate_sequences)
                for(auto &gene_b : it.second)
                    this->paralogues_best_hits.emplace_back(std::make_tuple(it.first, gene_b.first, gene_b.second));

            gene_id_pair.clear();
        }
    }

     bool check_constraint(int gene_id_a, int gene_id_b) {
        std::string sequence_a = this->sequences.operator[](gene_id_a);
        std::string sequence_b = this->sequences.operator[](gene_id_b);

         if(sequence_a.length() >= sequence_b.length()*2 || sequence_b.length() >= sequence_a.length()*2)
            return false;

        return true;
    }
};


#endif //PANDELOS_PLUSPLUS_PARALOGUES_H
