#ifndef PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H
#define PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "../../conf/conf.h"

class BidirectionalBestHits {
public:
    explicit BidirectionalBestHits(const std::unordered_map<int, std::unordered_map<int, double>>& best_hits, std::ofstream* log_stream) {
        if(log_stream->bad())
            throw std::runtime_error("the log stream's badbit error state flag is set");

        if(best_hits.empty())
            throw std::runtime_error("best hits map is empty");

        this->best_hits = &best_hits;
        this->log_stream = log_stream;
    }


    /*
     *  Date:
     *  A chiave di map_best_hits   (gene A)
     *  B valore di map_best_hits   (mappa C)
     *  D chiave di C               (gene D)
     *  E valore di C               (double)
     *
     *  Schema: map_best_hits< A, C<D, E> >
     *
     *  Si definisce un BIDIRECTIONAL BEST HITS, se esiste, una coppia di geni tale per cui:
     *
     *  per una chiave A, nella sua mappa C esiste un gene D che visto come una delle chiavi di map_best_hits,
     *  contiene nella sua mappa C una chiave uguale ad A.
    */
    void find_bidirectional_best_hits() {
        for (auto &map : *this->best_hits) {
            int id_gene_a = map.first;

            for(auto &submap : map.second) {
                int id_gene_b = submap.first;
                double jaccard = submap.second;

                //it: chiave uguale a gene_b
                auto it = this->best_hits->find(id_gene_b);

                if(it != this->best_hits->end()) {

                    //cerca il gene a come chiave della mappa di it (it->second)
                    auto it2 = it->second.find(id_gene_a);

                    if(it2 != it->second.end())
                        if(id_gene_a < id_gene_b)
                            this->bidirectional_best_hits.emplace_back(std::make_tuple(id_gene_a, id_gene_b, jaccard));
                }
            }
        }
    }

    const std::vector<std::tuple<int, int, double>>& get_bidirectional_best_hits() {
        if(this->bidirectional_best_hits.empty())
            throw std::runtime_error("bidirectional best hits vector is empty");

        return this->bidirectional_best_hits;
    }

private:
    std::ofstream* log_stream;
    const std::unordered_map<int, std::unordered_map<int, double>>* best_hits;
    std::vector<std::tuple<int, int, double>> bidirectional_best_hits;

};

#endif //PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H