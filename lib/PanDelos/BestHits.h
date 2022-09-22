#ifndef PANDELOS_PLUSPLUS_BESTHITS_H
#define PANDELOS_PLUSPLUS_BESTHITS_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "../../conf/conf.h"

class BestHits {
public:
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
