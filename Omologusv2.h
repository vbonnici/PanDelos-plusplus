#ifndef PANDELOS_PLUSPLUS_OMOLOGUSv2_H
#define PANDELOS_PLUSPLUS_OMOLOGUSv2_H

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
#include "include/global_options.h"


template <typename KeyType, typename LeftValue, typename RightValue>
std::map<KeyType, std::pair<LeftValue, RightValue>>
UnionMaps(const std::map<KeyType, LeftValue>& a1, const std::map<KeyType, RightValue>& a2) {
    std::map<KeyType, std::pair<LeftValue, RightValue>> result;
    typename std::map<KeyType, LeftValue>::const_iterator a1i  = a1.begin();
    typename std::map<KeyType, RightValue>::const_iterator a2i = a2.begin();

    while (a1i != a1.end() && a2i != a2.end()) {
        if (a1i->first < a2i->first) {
            result.insert(std::make_pair(a1i->first, std::make_pair(a1i->second, 0)));
            a1i++;
        } else if (a1i->first == a2i->first) {
            result.insert(std::make_pair(a1i->first, std::make_pair(a1i->second, a2i->second)));
            a1i++;
            a2i++;
        } else if (a1i->first > a2i->first) {
            result.insert(std::make_pair(a2i->first, std::make_pair(0, a2i->second)));
            a2i++;
        }
    }

    while (a1i != a1.end()) {
        result.insert(std::make_pair(a1i->first, std::make_pair(a1i->second, 0)));
        a1i++;
    }

    while (a2i != a2.end()) {
        result.insert(std::make_pair(a2i->first, std::make_pair(0, a2i->second)));
        a2i++;
    }

    return result;
}

class Omologusv2 {
public:
    explicit Omologusv2(const std::vector<std::string>&sequences,
                      const std::vector<std::pair<int, int>>&gene_pair_input,
                      const int sequences_type,
                      const int kmer_size,
                      std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size){

        this->log_stream = log_stream;
        this->gene_pair = &gene_pair_input;
        this->sequences_prefilter = &sequences;

        this->collect_sequences_id_from_gene_pair(&gene_pair_input);
    }

    /*
     * initializes a map for each sequence with the first k-mer and the counter at 0
     */
    void init_sequences_kmers() {

        for(auto &sequence_id : this->sequences_id) {

            std::map<std::bitset<kvalue>, int, bitset_comparer<kvalue>> temp;
            std::bitset<kvalue> bitset_temp;
            bitset_temp.set(false);
            temp.insert(std::make_pair(bitset_temp, 0));

            this->sequences_kmers.insert(std::make_pair(sequence_id, temp));
        }
    }

    void calculate_kmer_multiplicity() {
        std::string kmer;
        for(auto &i : this->sequences_kmers) {
            std::string sequence = this->sequences_prefilter->operator[](i.first);

            for(int window = 0; window < sequence.length(); ++window) {
                if(this->sequences_type == 1)
                    kmer = sequence.substr(window, this->kmer_size);
                else {
                    std::string aminoacid = sequence.substr(window, this->kmer_size); // /3

                    if(!aminoacid_is_valid(aminoacid))
                        continue;

                    kmer = Omologusv2::aminoacid_to_nucleotides(aminoacid);
                }

                if(kmer_is_valid(kmer)) {
                    std::bitset<kvalue> kmer_in_bit = kmer_to_bit(kmer);

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }
        }

        //*this->log_stream << "2 - kmer_multiplicity calcolate" << std::endl;
    }

    void calculate_best_hits(double jaccard_threshold = 0.0) {
        int id_gene_a;
        int id_gene_b;
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;

        //non inserisce duplicati
        for (auto &i: *this->gene_pair) {
            this->map_hits.insert(std::make_pair(i.first, std::unordered_map<int, double>()));
            this->map_hits.insert(std::make_pair(i.second, std::unordered_map<int, double>()));
        }

        if(!debug)
            omp_set_num_threads(omp_get_num_procs());

        #pragma omp parallel for private(jaccard_similarity, counter_min, counter_max, value_a, value_b, id_gene_a, id_gene_b)
        for (auto &i: *this->gene_pair) {
            id_gene_a = i.first;
            id_gene_b = i.second;

            counter_min = 0;
            counter_max = 0;

            auto kmer_key_a = this->sequences_kmers.find(id_gene_a);
            auto kmer_key_b = this->sequences_kmers.find(id_gene_b);

            std::map<std::string, int> temp_map1;

            auto map_a = kmer_key_a->second;

            for(auto &f : map_a)
                temp_map1.insert(std::make_pair(f.first.to_string(), f.second));

            std::map<std::string, int> temp_map2;

            auto map_b = kmer_key_b->second;

            for(auto &g : map_b)
                temp_map2.insert(std::make_pair(g.first.to_string(), g.second));

            auto result = UnionMaps<std::string, int, int>(temp_map1, temp_map2);

            for(auto &h : result) {

                auto y = h.second;

                value_a = y.first;
                value_b = y.second;

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
            //*this->log_stream << id_gene_a << " " << id_gene_b << " omologus jaccard similarity " << jaccard_similarity << std::endl;

            if (counter_max > 0 && jaccard_similarity >= jaccard_threshold) {
                auto result_a = this->map_hits.find(id_gene_a);

                #pragma omp critical
                result_a->second.insert(std::make_pair(id_gene_b, jaccard_similarity));


                auto result_b = this->map_hits.find(id_gene_b);

                #pragma omp critical
                result_b->second.insert(std::make_pair(id_gene_a, jaccard_similarity));
            }
        }

        //*this->log_stream << "2 - compute best hits " << std::endl;

    }

    std::vector<int>& get_sequences_id() {
        return this->sequences_id;
    }

    std::unordered_map<int, std::map<std::bitset<kvalue>, int, bitset_comparer<kvalue>>>& get_sequences_kmers() {
        return this->sequences_kmers;
    }

    std::unordered_map<int, std::unordered_map<int, double>> get_map_hits() {
        return this->map_hits;
    }

    std::vector<std::tuple<int, int, double>>& get_vector_tuple_bbh() {
        return this->vector_tuple_bbh;
    }

private:
    std::ofstream* log_stream;
    int kmer_size;
    const int sequences_type; //0 amino acids, 1 nucleotides
    std::vector<int> sequences_id;
    std::unordered_map<int, std::map<std::bitset<kvalue>, int, bitset_comparer<kvalue>>> sequences_kmers;   //map<sequence_id - map<kmer, contatore>>

    const std::vector<std::string>* sequences_prefilter;
    const std::vector<std::pair<int, int>>* gene_pair;
    std::unordered_map<int, std::unordered_map<int, double>> map_hits;
    std::vector<std::tuple<int, int, double>> vector_tuple_bbh;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        if(this->sequences_type == 0)
            return str.length() == this->kmer_size*3 && str.find_first_not_of("ACGT") == std::string::npos;
        else
            return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos;
    }

    [[nodiscard]] static bool aminoacid_is_valid(const std::string &str) {
        //return str.find_first_not_of("FLIMVSPTAY*HQNKDECWRG") == std::string::npos;
        return true;
    }

    static std::bitset<kvalue> kmer_to_bit(std::string& kmer) {
        std::bitset<kvalue> kmer_bit;
        kmer_bit.set(false);
        std::bitset<kvalue> A("00");
        std::bitset<kvalue> C("01");
        std::bitset<kvalue> G("10");
        std::bitset<kvalue> T("11");

        for(int a = 0; a < kmer.length(); ++a) {

            kmer_bit = kmer_bit << 2;

            if(reinterpret_cast<char>(kmer[a]) == 'A')
                kmer_bit = kmer_bit |= A;

            else if(reinterpret_cast<char>(kmer[a]) == 'C')
                kmer_bit = kmer_bit |= C;

            else if(reinterpret_cast<char>(kmer[a]) == 'G')
                kmer_bit = kmer_bit |= G;

            else if (reinterpret_cast<char>(kmer[a]) == 'T')
                kmer_bit = kmer_bit |= T;
        }

        return kmer_bit;
    }

    void collect_sequences_id_from_gene_pair(const std::vector<std::pair<int, int>>* gene_pair_input) {

        std::set<int> temp_sequences;

        for(auto &i : *gene_pair_input) {
            temp_sequences.insert(i.first);
            temp_sequences.insert(i.second);
        }

        this->sequences_id.reserve(temp_sequences.size());
        this->sequences_id.assign(temp_sequences.begin(), temp_sequences.end());
    }

    static std::string aminoacid_to_nucleotides(std::basic_string<char> aminoacid) {
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
#endif //PANDELOS_PLUSPLUS_OMOLOGUSv2_H
