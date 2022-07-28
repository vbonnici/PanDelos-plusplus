#ifndef PANDELOS_PLUSPLUS_PARALOG_H
#define PANDELOS_PLUSPLUS_PARALOG_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include "include/kvalue.h"
#include "include/global_options.h"
#include "Omologus.h"
#include "Omologusv2.h"
#include "PreFilter.h"

class Paralog {
public:
    explicit Paralog(std::vector<std::string> sequences_prefilter,
                     std::vector<std::vector<int>>& genome_sequencesid,
                     const int sequences_type,
                     int kmer_size,
                     std::vector<std::tuple<int, int, double>>& vector_tuple_bbh,
                     std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size) {

        this->log_stream = log_stream;
        this->sequences = std::move(sequences_prefilter);
        this->genome_sequencesid = genome_sequencesid;
        this->vector_tuple_bbh = &vector_tuple_bbh;

        this->genome_counter = this->genome_sequencesid.size();
    }


    void calculate_paralog() {
        this->genome_minimum_jaccard.reserve(this->genome_counter);
        this->max_genes_id.reserve(this->genome_counter);
        this->min_genes_id.reserve(this->genome_counter);

        for(int i = 0; i < this->genome_counter; ++i)
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
    int genome_counter;
    std::vector<std::vector<int>> genome_sequencesid;
    std::vector<double> genome_minimum_jaccard;
    std::vector<std::string> sequences;
    std::vector<std::tuple<int, int, double>>* vector_tuple_bbh;
    std::vector<std::tuple<int, int, double>> paralog_best_hits;
    std::vector<int> max_genes_id;
    std::vector<int> min_genes_id;


    void compute_minimum_jaccard() {
        for(int i = 0; i < this->genome_counter; ++i) {
            auto genes = this->genome_sequencesid.operator[](i);

            auto max = max_element(std::begin(genes), std::end(genes));
            auto min = min_element(std::begin(genes), std::end(genes));

            this->max_genes_id.operator[](i) = *max;
            this->min_genes_id.operator[](i) = *min;
        }

        for(int i = 0; i < this->genome_counter; ++i) {
            for(auto &tuple : *this->vector_tuple_bbh) {
                if(std::get<0>(tuple) >= this->min_genes_id.operator[](i) && std::get<0>(tuple) <= this->max_genes_id.operator[](i)) {
                    if(std::get<2>(tuple) < this->genome_minimum_jaccard.operator[](i))
                        this->genome_minimum_jaccard.operator[](i) = std::get<2>(tuple);
                }

                if(std::get<1>(tuple) >= this->min_genes_id.operator[](i) && std::get<1>(tuple) <= this->max_genes_id.operator[](i)) {
                    if(std::get<2>(tuple) < this->genome_minimum_jaccard.operator[](i))
                        this->genome_minimum_jaccard.operator[](i) = std::get<2>(tuple);
                }
            }
        }

        for(int i = 0; i < this->genome_counter; ++i) {
            *this->log_stream << "genome id " << i << " minimum jaccard " << this->genome_minimum_jaccard.operator[](i) << std::endl;
        }
    }

    void calculate_paralog_best_hits() {
        int max_gene_id;
        int min_gene_id;

        for(int i = 0; i < this->genome_counter; ++i) {
            max_gene_id = this->max_genes_id.operator[](i);
            min_gene_id = this->min_genes_id.operator[](i);

            std::vector<std::pair<int, int>> gene_id_pair;
            gene_id_pair.reserve(max_gene_id*max_gene_id);

            for(int a = min_gene_id; a <= max_gene_id; ++a)
                for(int b = a+1; b <= max_gene_id; ++b)
                    if(this->check_constraint(a, b))
                        gene_id_pair.emplace_back(std::make_pair(a, b));


            PreFilter prefilter = PreFilter(this->sequences, this->genome_sequencesid, this->sequences_type, this->log_stream);
            prefilter.init_sequences_kmers();
            prefilter.calculate_kmer_multiplicity();
            prefilter.calculate_best_hits(gene_id_pair, this->genome_minimum_jaccard.operator[](i));
            auto prefilter_best_hits = prefilter.get_best_hits();

            Omologusv2 omologus = Omologusv2(this->sequences, prefilter_best_hits, this->sequences_type, this->kmer_size, this->log_stream);

            omologus.init_sequences_kmers();
            omologus.calculate_kmer_multiplicity();
            omologus.calculate_best_hits(this->genome_minimum_jaccard.operator[](i));

            auto map_hits = omologus.get_map_hits();

            for(auto &it : map_hits)
                for(auto &gene_b : it.second)
                    this->paralog_best_hits.emplace_back(std::make_tuple(it.first, gene_b.first, gene_b.second));

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


#endif //PANDELOS_PLUSPLUS_PARALOG_H
