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
    explicit PreFiltering(const int& k, const std::vector<std::string>& sequences) : kmer_size(k){
        this->sequences = &sequences;
    }

    void populate_map_sequences() {
        for(auto &i: *this->sequences) {

            std::unordered_map<std::bitset<12>, int> temp;

            for (int j = 0; j < 1<<12 ; j++)
                temp.insert(std::make_pair(j, 0));                  //inizializza una mappa per ogni sequenza con tutte le possibili combinazioni di 12 bit

            this->map_sequences.insert(std::make_pair(i, temp)); //TODO: tempo di esecuzione molto lungo
        }
    }

    void calculate_kmer_frequency() {
        for(auto &i : *this->sequences) {
            for(int b = 0; b < i.length(); ++b)
                std::cout << i.substr(b, 6) << std::endl;

            //TODO:
            //piuttosto che calcolare tutte le combinazioni subito, inserisci solo 1 000.. e poi le altre
            //le inserisci o le aggiorni in questa funzione. Meglio usare array raw in vista della parallelizzazione

            //TODO:
            // - verificare che sia una sottostringa valida
            // - convertire in binario la sottostringa
            // - trovare la sequenza nella mappa
            // - aggiornare il contatore del k-mer
        }
    }

    std::unordered_map<std::string, std::unordered_map<std::bitset<12>, int>>& get_map_sequences() {
        return this->map_sequences;
    }

private:
    const int kmer_size; //TODO: per futuri utilizzi nucleotidico/amminoacidico - flag di riferimento
    const std::vector<std::string>* sequences;
    std::unordered_map<std::string, std::unordered_map<std::bitset<12>, int>> map_sequences; //<sequenza - <kmer, contatore>>

    bool string_is_valid(const std::string &str) {
        return str.find_first_not_of("ACGT") == std::string::npos;
    }
};

#endif //PANDELOS_PLUSPLUS_PREFILTERING_H
