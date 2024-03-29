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

/*
 * Class that allows the identification of paralogous genes. That is, the homology that occurs within the same genome
 */
class Paralogues {
public:

    /*
     * The constructor constructs the adequate data structures in order to proceed with the analysis to find the
     * paralogues genes
     *
     * @param[in] const std::vector<std::string>
     * @param[in] const std::vector<std::vector<int>>
     * @param[in] const std::vector<std::pair<int, int>>&
     * @param[in] const std::vector<std::tuple<int, int, double>>&
     * @param[in] const int
     * @param[in] const int
     * @param[in] std::ofstream*
     *
     * @throws std::runtime_error if the log stream's badbit error state flag is set
     * @throws std::runtime_error if the input sequences vector is empty
     * @throws std::runtime_error if the genome sequencesid vector is empty
     * @throws std::runtime_error if the genes id interval vector is empty
     * @throws std::runtime_error if the bidirectional best hits vector is empty
     */
    explicit Paralogues(const std::vector<std::string> input_sequences,
                        const std::vector<std::vector<int>> genome_sequencesid,
                        const std::vector<std::pair<int, int>>& genes_id_interval,
                        const std::vector<std::tuple<int, int, double>>& bidirectional_best_hits,
                        const int sequences_type,
                        const int kmer_size,
                        std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size) {

        if(log_stream->bad())
            throw std::runtime_error("the log stream's badbit error state flag is set");

        if(input_sequences.empty())
            throw std::runtime_error("input sequences vector is empty");

        if(genome_sequencesid.empty())
            throw std::runtime_error("genome sequencesid vector is empty");

        if(genes_id_interval.empty())
            throw std::runtime_error("genes id interval vector is empty");

        if(bidirectional_best_hits.empty())
            throw std::runtime_error("bidirectional best hits vector is empty");

        this->log_stream = log_stream;
        this->input_sequences = input_sequences;
        this->genome_sequencesid = genome_sequencesid;
        this->genes_id_interval = &genes_id_interval;
        this->bidirectional_best_hits = &bidirectional_best_hits;
        this->genome_counter = this->genome_sequencesid.size();
    }

    /*
     * This method initializes the minimum Jaccard similarity value, for each genome, to 1.0
     * (corresponding to the maximum possible value).
     *
     * He then proceeds to call the private methods which are in charge of calculating, for each genome,
     * the minimum value of actual Jaccard similarity and of identifying paralogous genes.
     */
    void find_paralogues() {
        for(int i = 0; i < this->genome_counter; ++i)
            this->genome_minimum_jaccard.push_back(1.0);

        this->calculate_minimum_jaccard();
        this->find_paralogues_best_hits();
    }

    /*
     * Getter that returns the vector of paralogous genes
     *
     * @param[out] const std::vector<std::tuple<int, int, double>>&
     */
    const std::vector<std::tuple<int, int, double>>& get_paralogues_best_hits() {
        return this->paralogues_best_hits;
    }

private:
    std::ofstream* log_stream;
    const int sequences_type;
    const int kmer_size;
    int genome_counter;
    std::vector<std::string> input_sequences;
    std::vector<std::vector<int>> genome_sequencesid;
    const std::vector<std::pair<int, int>>* genes_id_interval;
    const std::vector<std::tuple<int, int, double>>* bidirectional_best_hits;
    std::vector<double> genome_minimum_jaccard;
    std::vector<std::tuple<int, int, double>> paralogues_best_hits;


    /*
     * For each genome, for each occurrence of the bidirectional best hits genes, the minimum Jaccard similarity value is stored,
     * which will be used as the minimum threshold to be respected for the identification of paralogoues genes.
     *
     *
     */
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

    /*
     * To identify the paralogous genes we will use the classes already used to identify the orthologous genes, i.e.
     * the classes PreFilter.h and Homologues.h
     *
     * In this way it is possible, initially, to filter the input sequences and later it is possible to search for
     * paralogous genes considering higher thresholds, such as the kmer size value calculated at the beginning
     * by the software and the smallest Jaccard similarity value for the genome to consider.
     */
    void find_paralogues_best_hits() {
        int min_gene_id;
        int max_gene_id;

        for(int i = 0; i < this->genome_counter; ++i) {
            min_gene_id = std::get<0>(this->genes_id_interval->operator[](i));
            max_gene_id = std::get<1>(this->genes_id_interval->operator[](i));

            std::vector<std::pair<int, int>> gene_id_pair;

            //Vector containing all possible gene pairs for a given genome
            for(int a = min_gene_id; a <= max_gene_id; ++a)
                for(int b = a+1; b <= max_gene_id; ++b)
                    if(this->check_constraint(a, b))
                        gene_id_pair.emplace_back(std::make_pair(a, b));

            std::vector<std::pair<int, int>> pre_filter_candidate_sequences;
            try {
                PreFilter pre_filter = PreFilter(this->input_sequences, this->genome_sequencesid, this->sequences_type, this->log_stream);
                pre_filter.find_candidate_sequences(gene_id_pair, this->genome_minimum_jaccard.operator[](i));
                pre_filter_candidate_sequences = pre_filter.get_candidate_sequences();

            } catch(std::exception const& e) {
                *this->log_stream << "Exception: " << e.what() << std::endl;
                exit(11);
            }

            std::unordered_map<int, std::unordered_map<int, double>> homologues_candidate_sequences;
            if(!pre_filter_candidate_sequences.empty()) {
                try {
                    Homologues homologues = Homologues(this->input_sequences, pre_filter_candidate_sequences,
                                                       this->sequences_type, this->kmer_size, this->log_stream);
                    homologues.find_candidate_sequences(this->genome_minimum_jaccard.operator[](i));
                    homologues_candidate_sequences = homologues.get_candidate_sequences();

                } catch (std::exception const &e) {
                    *this->log_stream << "Exception: " << e.what() << std::endl;
                    exit(11);
                }

                for(auto &it : homologues_candidate_sequences)
                    for(auto &gene_b : it.second)
                        this->paralogues_best_hits.emplace_back(std::make_tuple(it.first, gene_b.first, gene_b.second));
            }

            gene_id_pair.clear();
        }
    }

    /*
     * Method that can be used to define constraints that the sequences must respect, even before carrying out the k-mer
     * comparison analysis
     *
     * @param[in] const int
     * @param[in] const int
     *
     * @param[out] bool
     */
     [[nodiscard]] bool check_constraint(const int gene_id_a, const int gene_id_b) const {
        const std::string sequence_a = this->input_sequences.operator[](gene_id_a);
        const std::string sequence_b = this->input_sequences.operator[](gene_id_b);

         if(sequence_a.length() >= sequence_b.length()*2 || sequence_b.length() >= sequence_a.length()*2)
            return false;

        return true;
    }
};


#endif //PANDELOS_PLUSPLUS_PARALOGUES_H