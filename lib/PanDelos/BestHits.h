#ifndef PANDELOS_PLUSPLUS_BESTHITS_H
#define PANDELOS_PLUSPLUS_BESTHITS_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "../../conf/conf.h"

/*
 * Class that allows to identify Best Hits relationships between sets of genes candidates to be homologous
 */
class BestHits {
public:

    /*
     * The constructor constructs the adequate data structures in order to proceed with the analysis to find the
     * best hits genes
     *
     * @param[in] const std::unordered_map<int, std::unordered_map<int, double>>&
     * @param[in] const std::vector<std::pair<int, int>>&
     * @param[in] std::ofstream*
     *
     * @throws std::runtime_error if the log_stream's badbit error state flag is set
     * @throws std::runtime_error if the input sequences map is empty
     * @throws std::runtime_error if the genome sequencesid vector is empty
     */
    explicit BestHits(const std::unordered_map<int, std::unordered_map<int, double>>& input_sequences,
                      const std::vector<std::pair<int, int>>& genes_id_interval,
                      std::ofstream* log_stream) {

        if(log_stream->bad())
            throw std::runtime_error("the log stream's badbit error state flag is set");

        if(input_sequences.empty())
            throw std::runtime_error("input sequences map is empty");

        if(genes_id_interval.empty())
            throw std::runtime_error("genome sequencesid vector is empty");

        this->log_stream = log_stream;
        this->input_sequences = &input_sequences;
        this->genes_id_interval = &genes_id_interval;
    }

    /*
     * This method identifies the best hits relationships, for each genome in this way.
     * For a given genome for which it is necessary to find the best hits genes, the range of identifiers of genes
     * to which it is associated is considered and for each of these, within the data structure containing similarity
     * values between genes, those with the highest similarity value.
     *
     * Note that for the same (maximum) similarity value, there can be multiple occurrences of gene pairs
     * that are best hits with this value.
     */
    void find_best_hits() {
        for (auto &it: *this->input_sequences) {
            if(!it.second.empty()) {
                this->best_hits.insert(std::make_pair(it.first, std::unordered_map<int, double>()));

                for(int genome_id = 0; genome_id < this->genes_id_interval->size(); ++genome_id) {
                    double max_jaccard = 0;
                    int max_gene_jaccard_id = -1;

                    for(auto &i: it.second) {
                        if(i.first >= std::get<0>(this->genes_id_interval->operator[](genome_id)) && i.first <= std::get<1>(this->genes_id_interval->operator[](genome_id))) {
                            if (i.second > max_jaccard) {
                                max_jaccard = i.second;
                                max_gene_jaccard_id = i.first;
                            }
                        }
                    }

                    if(max_gene_jaccard_id != -1) {
                        auto result = this->best_hits.find(it.first);
                        result->second.insert(std::make_pair(max_gene_jaccard_id, max_jaccard));
                    }
                }
            }
        }
    }

    /*
     * Getter which returns the map containing the identified best hits relationships and the corresponding similarity value
     *
     * @param[out] const std::unordered_map<int, std::unordered_map<int, double>>&
     * @throws std::runtime_error if the best_hits map is empty
     */
    const std::unordered_map<int, std::unordered_map<int, double>>& get_best_hits() {
        if(this->best_hits.empty())
            throw std::runtime_error("best_hits map is empty");

        return this->best_hits;
    }

private:
    std::ofstream* log_stream;
    const std::unordered_map<int, std::unordered_map<int, double>>* input_sequences;
    const std::vector<std::pair<int, int>>* genes_id_interval;
    std::unordered_map<int, std::unordered_map<int, double>> best_hits;

};

#endif //PANDELOS_PLUSPLUS_BESTHITS_H