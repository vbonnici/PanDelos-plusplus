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
    explicit PreFiltering(const int& k, const std::vector<std::string>& sequences) : kmer_size(k) {
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
                    //std::cout << kmer_in_bit << std::endl;

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }

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

    std::bitset<12> kmer_to_bit(std::string& kmer) {
        std::bitset<12> kmer_bit("000000000000");
        std::bitset<12> A("00");
        std::bitset<12> C("01");
        std::bitset<12> G("10");
        std::bitset<12> T("11");

        for(int a = 0; a < kmer.length(); ++a) {

            if(reinterpret_cast<char>(kmer[a]) == 'A')
                kmer_bit = kmer_bit |= A;

            else if(reinterpret_cast<char>(kmer[a]) == 'C')
                kmer_bit = kmer_bit |= C;

            else if(reinterpret_cast<char>(kmer[a]) == 'G')
                kmer_bit = kmer_bit |= G;

            else if (reinterpret_cast<char>(kmer[a]) == 'T')
                kmer_bit = kmer_bit |= T;

            kmer_bit = kmer_bit << 2;
        }

        return kmer_bit;
    }

    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>>& get_map_sequences() {
        return this->map_sequences;
    }

private:
    const int kmer_size; //TODO: per futuri utilizzi nucleotidico/amminoacidico - flag di riferimento
    const std::vector<std::string>* sequences;
    std::unordered_map<std::string, std::map<std::bitset<12>, int, bitset_comparer<12>>> map_sequences; //<sequenza - <kmer, contatore>>

    bool kmer_is_valid(const std::string &str) {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos; //TODO: aggiungere test di validità per amminoacidi
    }
};

#endif //PANDELOS_PLUSPLUS_PREFILTERING_H
