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
    explicit BidirectionalBestHits(std::unordered_map<int, std::unordered_map<int, double>> &map_best_hits, std::ofstream* log_stream) {
        this->map_best_hits = &map_best_hits;
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
    void calculate_bbh() {
        for (auto &map : *this->map_best_hits) {
            int id_gene_a = map.first;

            for(auto &submap : map.second) {
                int id_gene_b = submap.first;
                double jaccard = submap.second;

                //it: chiave uguale a gene_b
                auto it = this->map_best_hits->find(id_gene_b);

                if(it != this->map_best_hits->end()) {

                    //cerca il gene a come chiave della mappa di it (it->second)
                    auto it2 = it->second.find(id_gene_a);

                    if(it2 != it->second.end())
                        if(id_gene_a < id_gene_b)
                            this->vector_tuple_bbh.emplace_back(std::make_tuple(id_gene_a, id_gene_b, jaccard));
                }
            }
        }
    }

    std::vector<std::tuple<int, int, double>>& get_vector_tuple_bbh() {
        return this->vector_tuple_bbh;
    }

private:
    std::ofstream* log_stream;
    std::unordered_map<int, std::unordered_map<int, double>>* map_best_hits;
    std::vector<std::tuple<int, int, double>> vector_tuple_bbh;

};

#endif //PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H
