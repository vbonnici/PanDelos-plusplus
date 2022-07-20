#ifndef PANDELOS_PLUSPLUS_BESTHITS_H
#define PANDELOS_PLUSPLUS_BESTHITS_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "include/kvalue.h"
#include "include/global_options.h"

class BestHits {
public:
    explicit BestHits(std::unordered_map<int, std::unordered_map<int, double>> &map_hits, std::ofstream* log_stream) {
        this->map_hits = &map_hits;
        this->log_stream = log_stream;
    }

    void compute_best_hits() {
        for (auto &it: *this->map_hits) {
            if(!it.second.empty()) {
                this->map_best_hits.insert(std::make_pair(it.first, std::unordered_map<int, double>()));

                double max_jaccard = 0;

                for(auto &i: it.second) {
                    if(i.second > max_jaccard)
                        max_jaccard = i.second;
                }

                for(auto &i : it.second) {
                    if(i.second == max_jaccard) {
                        auto result = this->map_best_hits.find(it.first);

                        result->second.insert(std::make_pair(i.first, i.second));
                    }
                }
            }
        }
    }

    std::unordered_map<int, std::unordered_map<int, double>> get_map_best_hits() { //&
        return *this->map_hits; //this->map_best_hits
    }

    void print_map_best_hits() {
        for(auto &a : this->map_best_hits) {
            if(a.first == 0)
            for (auto &b: a.second)
                *this->log_stream << a.first << " " << b.first << " " << b.second << std::endl;
        }
    }

private:
    std::ofstream* log_stream;
    std::unordered_map<int, std::unordered_map<int, double>>* map_hits;
    std::unordered_map<int, std::unordered_map<int, double>> map_best_hits;

};
#endif //PANDELOS_PLUSPLUS_BESTHITS_H
