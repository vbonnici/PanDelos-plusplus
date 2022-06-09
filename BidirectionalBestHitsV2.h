#ifndef PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITSV2_H
#define PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITSV2_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <algorithm>

//aggiungere bitset comparer

//typedef std::unordered_map<int, int>::value_type unmap_string_string_value_type;

//aggiungere StartWith

class BidirectionalBestHitsV2 {
public:
    explicit BidirectionalBestHitsV2(const std::vector<std::string>& sequences,
                                   std::vector<std::pair<int, int>>& best_hits,
                                   std::vector<std::vector<int>>& genome_sequenceid,
                                   const int flag) : flag(flag), jaccard_threshold(0.8) {
        this->compute_alphabet(&sequences);
        this->kmer_size = 9; //TODO: scelta di k dinamico
        this->collect_sequences_from_best_hits(&sequences, &best_hits);
        std::cout << "costruttore bbh" << std::endl;
        this->best_hits_prefilter = &best_hits;
        this->sequences_prefilter = &sequences;
        this->genome_sequenceid = &genome_sequenceid;
    }

    /*
     * initializes a map for each sequence with the first k-mer and the counter at 0
     */
    void init_sequences_kmers() {

        for(auto &sequence: this->sequences) {

            std::map<std::bitset<18>, int, bitset_comparer<18>> temp;

            temp.insert(std::make_pair("00000000000000000000", 0));

            this->sequences_kmers.insert(std::make_pair(sequence, temp));
        }
    }

    void calculate_kmer_multiplicity() {
        int counter = 0;
        std::string kmer;
        for(auto &i : this->sequences_kmers) {
            std::string sequence = i.first;

            for(int window = 0; window < sequence.length(); ++window) {
                if(this->flag == 1)
                    kmer = sequence.substr(window, this->kmer_size);
                else
                    kmer = BidirectionalBestHitsV2::aminoacid_to_nucleotides(sequence.substr(window, 3));

                //std::cout << kmer << std::endl;
                if(kmer_is_valid(kmer)) {
                    //std::cout << kmer << std::endl;
                    std::bitset<18> kmer_in_bit = kmer_to_bit(kmer);
                    //std::cout << kmer_in_bit << std::endl;

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }
            //++counter;
            //std::cout << sequence[0] << counter << std::endl;
        }

        std::cout << "kmer_multiplicity calcolate" << std::endl;
    }

    /*void calculate_best_hits() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;

        int index;

        for(auto &a : this->sequences_kmers) {
            for(auto &b : this->sequences_kmers) {
                sequence_a = a.first;
                sequence_b = b.first;

                counter_min = 0;
                counter_max = 0;

                if(sequence_a != sequence_b) { //TODO: aggiungere vincoli secondo calcolo bh
                    for (int j = 0; j < 1<<18 ; j++) {
                        std::bitset<18> kmer(j);
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

                        counter_min += std::min(value_a, value_b);
                        counter_max += std::max(value_a, value_b);
                    }


                    //std::cout << "min " << counter_min << " max " << counter_max << std::endl;

                    jaccard_similarity = (double) counter_min / counter_max;
                    //std::cout << "jaccard similarity " << jaccard_similarity << std::endl;

                    if(std::isfinite(jaccard_similarity) && jaccard_similarity > this->jaccard_threshold) {
                        std::unordered_map<std::string, double> temp;
                        temp.insert(std::make_pair(sequence_b, jaccard_similarity));

                        this->map_best_hits.insert(std::make_pair(sequence_a, temp));
                    }
                }
            }
        }
        std::cout << "seconda parte best hits calcolati " << std::endl;

    }*/

    void calculate_best_hits_v2() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;

        int index;

        for(index = 0; index < this->genome_sequenceid->size(); ++index) {
            for(int i = index + 1; i < this->genome_sequenceid->size(); ++i) {

                        std::vector<int> genome_a = this->genome_sequenceid->operator[](index);
                        std::vector<int> genome_b = this->genome_sequenceid->operator[](i);

                        for (auto &a: genome_a)
                            for (auto &b: genome_b) {

                                /*auto it = std::find_if(this->best_hits_prefilter->begin(),
                                                               this->best_hits_prefilter->end(), FindPair(a, b));*/

                                auto it = std::find_if( this->best_hits_prefilter->begin(), this->best_hits_prefilter->end(),
                                                        [&a, &b](const std::pair<int, int>& element){ return element.first == a && element.second == b;} );

                                    if (it != this->best_hits_prefilter->end()) {
                                        //std::cout << a << " " << b << std::endl;

                                        sequence_a = this->sequences_prefilter->operator[](a);
                                        sequence_b = this->sequences_prefilter->operator[](b);

                                        counter_min = 0;
                                        counter_max = 0;

                                        auto kmer_key_a = this->sequences_kmers.find(sequence_a);
                                        auto kmer_key_b = this->sequences_kmers.find(sequence_b);

                                        for (int j = 0; j < 1 << 18; j++) {
                                            std::bitset<18> kmer(j);
                                            //std::cout << kmer << std::endl;

                                            auto result_a = kmer_key_a->second.find(kmer);
                                            auto result_b = kmer_key_b->second.find(kmer);

                                            value_a = 0;
                                            value_b = 0;

                                            if (result_a != kmer_key_a->second.end()) {
                                                if (result_a->first == kmer) {
                                                    value_a = result_a->second;
                                                    //std::cout << "j " << j << " result_a->second " << result_a->second << std::endl;
                                                }
                                            }

                                            if (result_b != kmer_key_b->second.end()) {
                                                if (result_b->first == kmer) {
                                                    value_b = result_b->second;
                                                    //std::cout << "j " << j << " result_b->second " << result_b->second << " value_b " << value_b << std::endl;
                                                }
                                            }

                                            counter_min += std::min(value_a, value_b);
                                            counter_max += std::max(value_a, value_b);
                                        }

                                        jaccard_similarity = (double) counter_min / counter_max;
                                        //std::cout << "prefilter jaccard similarity " << jaccard_similarity << std::endl;

                                        if (counter_max > 0 && jaccard_similarity > this->jaccard_threshold) {
                                            this->vector_best_hits.emplace_back(std::make_pair(a, b));
                                            //std::cout << a << " best hit con " << b << std::endl;
                                            std::unordered_map<int, double> temp;
                                            temp.insert(std::make_pair(b, jaccard_similarity));

                                            this->new_map_best_hits.insert(std::make_pair(a, temp));
                                        }
                                    }
                            }

            }
        }

        std::cout << "seconda parte best hits calcolati v2 " << std::endl;
    }

    void calculate_bbh() {
        for(auto coppia : this->vector_best_hits) {
            int a = coppia.first;
            int b = coppia.second;
            auto it = std::find_if( this->vector_best_hits.begin(), this->vector_best_hits.end(),
                                    [&a, &b](const std::pair<int, int>& element){ return element.first == b && element.second == a;} );

            if(it != this->vector_best_hits.end())
                std::cout << it->first << " " << it->second << std::endl;
        }
    }

    /*
     * Per ogni occorrenza della mappa, se il gene è diverso, confrontiamo 2 occorrenze alla volta (seconda parte mappa)
     * Cerchiamo ogni occorrenza di questa seconda parte del PRIMO GENE nella seconda parte della mappa del SECONDO GENE
     * Poi lo facciamo al contrario. Se c'è corrispondenza aggiungiamo questo nuovo bbh alla lista
     *
     */
    void calculate_bidirectional_best_hits() {

        std::unordered_map<int, int> map_bidirectional_best_hits_internal;

        for(auto &gene_a: this->new_map_best_hits) {
            auto gene_a_best_hits = gene_a.second;

            for(auto &gene_b: this->new_map_best_hits) {
                auto gene_b_best_hits = gene_b.second;

                if(gene_a.first != gene_b.first) {

                    auto c = gene_b_best_hits.find(gene_a.first);
                    auto d = gene_a_best_hits.find(gene_b.first);

                    /// if true, gene_a is a best hit of gene_b and vicevers
                    if(c != gene_b_best_hits.end() && d != gene_a_best_hits.end()) {

                        ///check if a key with gene c already exists
                        auto temp_a = map_bidirectional_best_hits_internal.find(c->first);

                        ///check if a value with gene d already exists
                        auto temp_b = std::find_if(map_bidirectional_best_hits_internal.begin(),
                                                   map_bidirectional_best_hits_internal.end(),
                                                   [&d](const std::unordered_map<int, int>::value_type& vt)
                                                    {
                                                        return vt.second == d->first;
                                                    });
                        ///if true, it means that a key and a value already exist for genes c and d
                        if(temp_a != map_bidirectional_best_hits_internal.end() && temp_b != map_bidirectional_best_hits_internal.end())
                            continue;

                        ///perform the same procedure but perform a key-value search with the genes reversed

                        temp_a = map_bidirectional_best_hits_internal.find(d->first);

                        temp_b = std::find_if(map_bidirectional_best_hits_internal.begin(),
                                              map_bidirectional_best_hits_internal.end(),
                                                    [&c](const std::unordered_map<int, int>::value_type& vt)
                                                    {
                                                        return vt.second == c->first;
                                                    });

                        if(temp_a != map_bidirectional_best_hits_internal.end() && temp_b != map_bidirectional_best_hits_internal.end())
                            continue;

                        ///arrived here means that there is no risk of inserting a duplicate key-value or key-value
                        map_bidirectional_best_hits_internal.insert(std::make_pair(c->first, d->first));

                        std::unordered_map<unsigned int, double> temp;
                        temp.insert(std::make_pair(d->first, d->second));

                        std::cout << c->first << " " << d->first << std::endl;

                        this->map_bidirectional_best_hits.insert(std::make_pair(c->first, temp));
                    }
                }
            }
        }

        std::cout << "bbh calcolati" << std::endl;
    }

    /*
     * For each gene sequence, it inserts each character of which it is composed into an unordered set
     * (if not yet present)
     */
    void compute_alphabet(const std::vector<std::string>* sequences_input) {
        for (auto &i: *sequences_input) {
            for(char a : i)
                this->alphabet.insert(a);
        }

        std::cout << "alfabeto calcolato" << std::endl;
    }

    std::vector<std::string>& get_sequences() {
        return this->sequences;
    }

    std::unordered_map<std::string, std::map<std::bitset<18>, int, bitset_comparer<18>>>& get_sequences_kmers() {
        return this->sequences_kmers;
    }

    std::unordered_map<std::string, std::unordered_map<std::string, double>>& get_map_best_hits() {
        return this->map_best_hits;
    }

    std::unordered_map<int, std::unordered_map<int, double>>& get_new_map_best_hits() {
        return this->new_map_best_hits;
    }

    std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>>& get_map_bidirectional_best_hits() {
        return this->map_bidirectional_best_hits;
    }

    /*
     * Getter which returns the unsorted set containing the alphabet obtained from all genes in the .faa file
     *
     * @param[out] std::unordered_set<char>&
     */
    std::unordered_set<char>& get_alphabet() {
        return this->alphabet;
    }

    std::vector<std::pair<int, int>> get_vector_best_hits() {
        return this->vector_best_hits;
    }

private:
    int kmer_size;
    const int flag; //0 amino acids, 1 nucleotides
    const double jaccard_threshold;
    std::vector<std::string> sequences;
    const std::vector<std::string>* sequences_prefilter;
    std::unordered_map<std::string, std::map<std::bitset<18>, int, bitset_comparer<18>>> sequences_kmers; //map<sequenza - map<kmer, contatore>>
    std::vector<std::pair<int, int>>* best_hits_prefilter;
    std::vector<std::vector<int>>* genome_sequenceid;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> map_best_hits; //map<sequence_name1, map<sequence_name2, jaccard>>
    std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> map_bidirectional_best_hits; //map<sequence_name1, map<sequence_name2, jaccard>>
    std::unordered_set<char> alphabet;
    std::vector<std::pair<int, int>> vector_best_hits;
    std::unordered_map<int, std::unordered_map<int, double>> new_map_best_hits;
    std::tuple<int, int, double> tuple_bidirectional_best_hits;

    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos; //TODO: aggiungere test di validità per amminoacidi
    }

    static std::bitset<18> kmer_to_bit(std::string& kmer) {
        std::bitset<18> kmer_bit("00000000000000000000");
        std::bitset<18> A("00");
        std::bitset<18> C("01");
        std::bitset<18> G("10");
        std::bitset<18> T("11");

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

    static bool check_constraint(std::string& sequence_a, std::string& sequence_b) {
        if(sequence_a == sequence_b)
            return false;

        if(sequence_a.length() >= sequence_b.length()*2 || sequence_b.length() >= sequence_a.length()*2)
            return false;

        return true;
    }

    void collect_sequences_from_best_hits(const std::vector<std::string>* sequences_input,
                                          std::vector<std::pair<int, int>>* best_hits_input) {

        std::set<std::string> temp_sequences;

        for(auto &i : *best_hits_input) {
            std::string sequence_a = sequences_input->operator[](i.first);
            std::string sequence_b = sequences_input->operator[](i.second);
            temp_sequences.insert(sequence_a);
            temp_sequences.insert(sequence_b);
        }

        this->sequences.reserve(temp_sequences.size());
        this->sequences.assign(temp_sequences.begin(), temp_sequences.end());

        std::cout << "sequences collected" << std::endl;
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

#endif //PANDELOS_PLUSPLUS_BIDIRECTIONALBESTHITSV2_H
