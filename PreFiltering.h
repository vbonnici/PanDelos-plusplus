#ifndef PANDELOS_PLUSPLUS_PREFILTERING_H
#define PANDELOS_PLUSPLUS_PREFILTERING_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <bitset>
#include <unordered_map>


template<size_t sz> struct bitset_comparer { //serve solo per le map ordinate
    bool operator() (const std::bitset<sz> &b1, const std::bitset<sz> &b2) const {
        return b1.to_ulong() < b2.to_ulong();
    }
};

class PreFiltering {
public:
    explicit PreFiltering(const int& k, const std::vector<std::string>& sequences) : kmer_size(k), soglia_jaccard(0.8) {
        this->sequences = &sequences;
    }

    void populate_map_sequences() {
        for(auto &i: *this->sequences) {

            std::map<std::bitset<12>, int, bitset_comparer<12>> temp;

            temp.insert(std::make_pair("000000000000", 0));                  //inizializza una mappa per ogni sequenza con il primo k-mer e il contatore a 0

            this->map_sequences.insert(std::make_pair(i, temp));
        }
    }

    void calculate_kmer_frequency() {
        for(auto &i : this->map_sequences) {
            std::string sequence = i.first;
            for(int b = 0; b < sequence.length(); ++b) {
                std::string kmer = sequence.substr(b, this->kmer_size);
                if(kmer_is_valid(kmer)) {
                    //std::cout << kmer << std::endl;
                    std::bitset<12> kmer_in_bit = kmer_to_bit(kmer);
                    std::cout << kmer_in_bit << std::endl;

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }
        }
    }

    void calculate_bh() {
        for(auto &a : this->map_sequences) {
            for(auto &b : this->map_sequences) {
                std::vector<std::string> sequence_a = split_string(a.first, ',');
                std::vector<std::string> sequence_b = split_string(b.first, ',');

                int value_min[4095];
                int value_max[4095];
                int index = 0;

                if(sequence_a[0] != sequence_b[0]) {
                    for (int j = 0; j < 1<<12 ; j++) {
                        std::bitset<12> kmer(j);

                        auto result_a = a.second.find(kmer);
                        auto result_b = b.second.find(kmer);

                        int value_a = 0;
                        int value_b = 0;

                        if(result_a != a.second.end())
                            value_a = result_a->second;

                        if(result_b != b.second.end())
                            value_b = result_b->second;

                        if(value_a < value_b) {
                            value_min[index] = value_a;
                            ++index;
                        }
                        else {
                            value_max[index] = value_b;
                            ++index;
                        }

                        double jaccard_similarity = sum_array(value_min, 4095) / sum_array(value_max, 4095);

                        if(jaccard_similarity > this->soglia_jaccard) {
                            std::unordered_map<std::string, double> temp;
                            temp.insert(std::make_pair(sequence_b[0], jaccard_similarity));

                            this->map_genes_jaccard.insert(std::make_pair(sequence_a[0], temp));
                        }
                    }
                }
            }
        }
    }

    /*
     * Per ogni occorrenza della mappa, se il gene è diverso, confrontiamo 2 occorrennze alla volta (seconda parte mappa)
     * Cerchiamo ogni occorrenza di questa seconda parte del PRIMO GENE nella seconda parte della mappa del SECONDO GENE
     * Poi lo facciamo al contrario. Se c'è corrispondenza aggiungiamo questo nuovo bbh alla lista
     *
     */
    void calculate_bbh() {
        for(auto &a: this->map_genes_jaccard) {
            auto second_genes_map_a = a.second;
            for(auto &b: this->map_genes_jaccard) {
                auto second_genes_map_b = b.second;

                if(a.first != b.first) {
                    for (auto &c: second_genes_map_a) {
                        auto result_a = second_genes_map_b.find(c.first);
                        if (result_a != second_genes_map_b.end()) {
                            for (auto &d: second_genes_map_b) {
                                auto result_b = second_genes_map_a.find(d.first);

                                if (result_b != second_genes_map_a.end()) {
                                    this->genes_bbh.insert(std::make_pair(result_a->first, result_b->first));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::bitset<12> kmer_to_bit(std::string& kmer) {
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

    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>>& get_map_sequences() {
        return this->map_sequences;
    }

private:
    const int kmer_size; //TODO: per futuri utilizzi nucleotidico/amminoacidico - flag di riferimento
    const double soglia_jaccard;
    const std::vector<std::string>* sequences;
    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>> map_sequences; //<sequenza - <kmer, contatore>>
    std::unordered_map<std::string, std::unordered_map<std::string, double>> map_genes_jaccard; //<sequence_name1, <sequence_name2, jaccard>>
    std::unordered_map<std::string, std::string> genes_bbh; //<sequence_name1, <sequence_name2>

    bool kmer_is_valid(const std::string &str) {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos; //TODO: aggiungere test di validità per amminoacidi
    }

    std::vector<std::string> split_string(std::string const &str, const char delim) {
        std::vector<std::string> out;
        size_t start;
        size_t end = 0;

        while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
            end = str.find(delim, start);
            out.push_back(str.substr(start, end - start));
        }

        return out;
    }

    static int sum_array(const int arr[], int n)
    {
        int sum = 0;

        for (int i = 0; i < n; i++)
            sum += arr[i];

        return sum;
    }
};

#endif //PANDELOS_PLUSPLUS_PREFILTERING_H
