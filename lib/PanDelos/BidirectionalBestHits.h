#ifndef PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H
#define PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITS_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "../../conf/conf.h"

/*
 * Class that allows to identify Bidirectional Best Hits relationships between sets of genes which are best hits.
 */
class BidirectionalBestHits {
public:

    /*
     * The constructor constructs the adequate data structures in order to proceed with the analysis to find the
     * bidirectional best hits genes
     *
     * @param[in] const std::unordered_map<int, std::unordered_map<int, double>>&
     * @param[in] std::ofstream*
     *
     * @throws std::runtime_error if the log_stream's badbit error state flag is set
     * @throws std::runtime_error if the best hits map is empty
     */
    explicit BidirectionalBestHits(const std::unordered_map<int, std::unordered_map<int, double>>& best_hits, std::ofstream* log_stream) {
        if(log_stream->bad())
            throw std::runtime_error("the log stream's badbit error state flag is set");

        if(best_hits.empty())
            throw std::runtime_error("best hits map is empty");

        this->best_hits = &best_hits;
        this->log_stream = log_stream;
    }


    /*
     *  This definition of Bidirectional Best Hits is essential and self-explanatory to understand how the method works.
     *
     *  Having:
     *
     *  A key of map_best_hits   (gene A)
     *  B value of map_best_hits (map C)
     *  D key of C               (gene D)
     *  E value of C             (double)
     *
     *  Scheme: map_best_hits <A, C <D, E>>
     *
     *  A BIDIRECTIONAL BEST HITS is defined, a pair of genes such that:
     *
     *  for a key A, in its map C there is a gene D which, seen as one of the keys of map_best_hits,
     *  contains in its map C a key equal to A.
    */
    void find_bidirectional_best_hits() {
        for (auto &map : *this->best_hits) {
            int id_gene_a = map.first;

            for(auto &submap : map.second) {
                int id_gene_b = submap.first;
                double jaccard = submap.second;

                //it: key equal to gene_b
                auto it = this->best_hits->find(id_gene_b);

                if(it != this->best_hits->end()) {

                    //search for gene a as key to it map (it-> second)
                    auto it2 = it->second.find(id_gene_a);

                    if(it2 != it->second.end())
                        if(id_gene_a < id_gene_b)
                            this->bidirectional_best_hits.emplace_back(std::make_tuple(id_gene_a, id_gene_b, jaccard));
                }
            }
        }
    }

    /*
     * Getter which returns the vector containing the identified bidirectional best hits relationships and the corresponding similarity value
     *
     * @param[out] const std::vector<std::tuple<int, int, double>>&
     * @throws std::runtime_error if the bidirectional best hits vector is empty
     */
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