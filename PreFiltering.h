#ifndef PANDELOS_PLUSPLUS_PREFILTERING_H
#define PANDELOS_PLUSPLUS_PREFILTERING_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>


template<size_t sz> struct bitset_comparer {
    bool operator() (const std::bitset<sz> &b1, const std::bitset<sz> &b2) const {
        return b1.to_ulong() < b2.to_ulong();
    }
};

typedef std::unordered_map<std::string, std::string>::value_type unmap_string_string_value_type;


class PreFiltering {
public:
    explicit PreFiltering(const int& k, const std::vector<std::string>& sequences) : kmer_size(k), jaccard_threshold(0.8) {
        this->sequences = &sequences;
    }

    /*
     * initializes a map for each sequence with the first k-mer and the counter at 0
     */
    void init_map_sequences_kmers() {
        for(auto &sequence: *this->sequences) {

            std::map<std::bitset<12>, int, bitset_comparer<12>> temp;

            temp.insert(std::make_pair("000000000000", 0));

            this->map_sequences_kmers.insert(std::make_pair(sequence, temp));
        }
    }

    void calculate_kmer_frequency() {
        for(auto &i : this->map_sequences_kmers) {
            std::vector<std::string> sequence = split_string(i.first, ',');

            for(int window = 0; window < sequence[1].length(); ++window) {
                std::string kmer = sequence[1].substr(window, this->kmer_size);

                if(kmer_is_valid(kmer)) {
                    //std::cout << kmer << std::endl;
                    std::bitset<12> kmer_in_bit = kmer_to_bit(kmer);
                    //std::cout << kmer_in_bit << std::endl;

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }
        }
    }

    //TODO: cancellare stampe di debug
    void calculate_best_hits() {
        int value_a;
        int value_b;
        int counter_min;
        int counter_max;
        double jaccard_similarity;

        for(auto &a : this->map_sequences_kmers) {
            for(auto &b : this->map_sequences_kmers) {
                std::vector<std::string> sequence_a = split_string(a.first, ',');
                std::vector<std::string> sequence_b = split_string(b.first, ',');

                counter_min = 0;
                counter_max = 0;

                if(sequence_a[0] != sequence_b[0]) {
                    for (int j = 0; j < 1<<12 ; j++) {
                        std::bitset<12> kmer(j);
                        //std::cout << kmer << std::endl;

                        auto result_a = a.second.find(kmer);
                        auto result_b = b.second.find(kmer);

                        value_a = 0;
                        value_b = 0;

                        if (result_a != a.second.end()) {
                            if(result_a->first == kmer) {
                                value_a = result_a->second;
                                //std::cout << "j " << j << " result_a->second " << result_a->second << std::endl;
                            }
                        }

                        if (result_b != b.second.end()) {
                            if(result_b->first == kmer) {
                                value_b = result_b->second;
                                //std::cout << "j " << j << " result_b->second " << result_b->second << " value_b " << value_b << std::endl;
                            }
                        }


                        if (value_a < value_b || value_a == 0) {
                            if(value_a == 1 && value_b == 2)
                                std::cout  << " value_a " << value_a << " value_b " << value_b<< std::endl;

                            counter_min += value_a;

                        } else {
                            if(value_a == 1 && value_b == 1)
                                std::cout  << " value_a " << value_a << " value_b " << value_b<< std::endl;

                            counter_max += value_a;
                        }
                    }


                    std::cout << "min " << counter_min << " max " << counter_max << std::endl;

                    jaccard_similarity = (double) counter_min / counter_max;
                    std::cout << "jaccard similarity " << jaccard_similarity << std::endl;

                    if(isfinite(jaccard_similarity) && jaccard_similarity > this->jaccard_threshold) {
                        std::unordered_map<std::string, double> temp;
                        temp.insert(std::make_pair(sequence_b[0], jaccard_similarity));

                        this->map_best_hits.insert(std::make_pair(sequence_a[0], temp));
                    }
                }
            }
        }
    }

    /*
     * Per ogni occorrenza della mappa, se il gene è diverso, confrontiamo 2 occorrenze alla volta (seconda parte mappa)
     * Cerchiamo ogni occorrenza di questa seconda parte del PRIMO GENE nella seconda parte della mappa del SECONDO GENE
     * Poi lo facciamo al contrario. Se c'è corrispondenza aggiungiamo questo nuovo bbh alla lista
     *
     */
    void calculate_bidirectional_best_hits() {
        for(auto &gene_a: this->map_best_hits) {
            auto gene_a_best_hits = gene_a.second;

            for(auto &gene_b: this->map_best_hits) {
                auto gene_b_best_hits = gene_b.second;

                if(gene_a.first != gene_b.first) {

                    auto c = gene_b_best_hits.find(gene_a.first);
                    auto d = gene_a_best_hits.find(gene_b.first);

                    /// if true, gene_a is a best hit of gene_b and vicevers
                    if(c != gene_b_best_hits.end() && d != gene_a_best_hits.end()) {

                        ///check if a key with gene c already exists
                        auto temp_a = this->map_bidirectional_best_hits.find(c->first);

                        ///check if a value with gene d already exists
                        auto temp_b = std::find_if(this->map_bidirectional_best_hits.begin(),
                                                   this->map_bidirectional_best_hits.end(),
                                                   [&d](const unmap_string_string_value_type& vt)
                                                    {
                                                        return vt.second == d->first;
                                                    });
                        ///if true, it means that a key and a value already exist for genes c and d
                        if(temp_a != this->map_bidirectional_best_hits.end() && temp_b != this->map_bidirectional_best_hits.end())
                            continue;

                        ///perform the same procedure but perform a key-value search with the genes reversed

                        temp_a = this->map_bidirectional_best_hits.find(d->first);

                        temp_b = std::find_if(this->map_bidirectional_best_hits.begin(),
                                              this->map_bidirectional_best_hits.end(),
                                                    [&c](const unmap_string_string_value_type& vt)
                                                    {
                                                        return vt.second == c->first;
                                                    });

                        if(temp_a != this->map_bidirectional_best_hits.end() && temp_b != this->map_bidirectional_best_hits.end())
                            continue;

                        ///arrived here means that there is no risk of inserting a duplicate key-value or key-value
                        this->map_bidirectional_best_hits.insert(std::make_pair(c->first, d->first));
                    }
                }
            }
        }
    }

    static std::bitset<12> kmer_to_bit(std::string& kmer) {
        std::bitset<12> kmer_bit("000000000000");
        std::bitset<12> A("00");
        std::bitset<12> C("01");
        std::bitset<12> G("10");
        std::bitset<12> T("11");

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

    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>>& get_map_sequences_kmers() {
        return this->map_sequences_kmers;
    }

    std::unordered_map<std::string, std::unordered_map<std::string, double>>& get_map_best_hits() {
        return this->map_best_hits;
    }

    std::unordered_map<std::string, std::string>& get_map_bidirectional_best_hits() {
        return this->map_bidirectional_best_hits;
    }

private:
    const int kmer_size; //TODO: per futuri utilizzi nucleotidico/amminoacidico - flag di riferimento
    const double jaccard_threshold;
    const std::vector<std::string>* sequences; //TODO: ricorda di cambiare nome...
    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>> map_sequences_kmers; //<sequenza - <kmer, contatore>>
    std::unordered_map<std::string, std::unordered_map<std::string, double>> map_best_hits; //<sequence_name1, <sequence_name2, jaccard>>
    std::unordered_map<std::string, std::string> map_bidirectional_best_hits; //<sequence_name1, <sequence_name2>

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos; //TODO: aggiungere test di validità per amminoacidi
    }

    static std::vector<std::string> split_string(std::string const &str, const char delim) {
        std::vector<std::string> out;
        size_t start;
        size_t end = 0;

        while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
            end = str.find(delim, start);
            out.push_back(str.substr(start, end - start));
        }

        return out;
    }

};

#endif //PANDELOS_PLUSPLUS_PREFILTERING_H
