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

class BidirectionalBestHitsV2 {
public:
    explicit BidirectionalBestHitsV2(const std::vector<std::string>& sequences,
                                   std::vector<std::pair<int, int>>& best_hits,
                                   std::vector<std::vector<int>>& genome_sequenceid,
                                   const int flag) : flag(flag) {
        this->compute_alphabet(&sequences);
        this->kmer_size = 9; //TODO: scelta di k dinamico
        this->collect_sequences_from_best_hits(&sequences, &best_hits);
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
        std::string kmer;
        for(auto &i : this->sequences_kmers) {
            std::string sequence = i.first;

            for(int window = 0; window < sequence.length(); ++window) {
                if(this->flag == 1)
                    kmer = sequence.substr(window, this->kmer_size);
                else
                    kmer = BidirectionalBestHitsV2::aminoacid_to_nucleotides(sequence.substr(window, 3));

                if(kmer_is_valid(kmer)) {
                    std::bitset<18> kmer_in_bit = kmer_to_bit(kmer);

                    auto result = i.second.find(kmer_in_bit);

                    if(result == i.second.end())
                        i.second.insert(std::make_pair(kmer_in_bit, 1));
                    else
                        result->second += 1;
                }
            }
        }

        std::cout << "2 - kmer_multiplicity calcolate" << std::endl;
    }

    void calculate_best_hits() {
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;
        std::string sequence_a;
        std::string sequence_b;

        for (auto &i: *this->best_hits_prefilter) {
            this->map_hits.insert(std::make_pair(i.first, std::unordered_map<int, double>()));
            this->map_hits.insert(std::make_pair(i.second, std::unordered_map<int, double>()));
        }

        for (auto &i: *this->best_hits_prefilter) {
            int id_gene_a = i.first;
            int id_gene_b = i.second;

            sequence_a = this->sequences_prefilter->operator[](id_gene_a);
            sequence_b = this->sequences_prefilter->operator[](id_gene_b);

            counter_min = 0;
            counter_max = 0;

            auto kmer_key_a = this->sequences_kmers.find(sequence_a);
            auto kmer_key_b = this->sequences_kmers.find(sequence_b);

            for (int j = 0; j < 1 << 18; j++) {
                std::bitset<18> kmer(j);

                auto result_a = kmer_key_a->second.find(kmer);
                auto result_b = kmer_key_b->second.find(kmer);

                value_a = 0;
                value_b = 0;

                if (result_a != kmer_key_a->second.end())
                    if (result_a->first == kmer)
                        value_a = result_a->second;


                if (result_b != kmer_key_b->second.end())
                    if (result_b->first == kmer)
                        value_b = result_b->second;

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
            //std::cout << "2 - jaccard similarity " << jaccard_similarity << std::endl;

            if (counter_max > 0) {
                this->map_hits[id_gene_a].insert(std::make_pair(id_gene_b, jaccard_similarity));
                this->map_hits[id_gene_b].insert(std::make_pair(id_gene_a, jaccard_similarity));
            }
        }

        this->compute_best_hits();

        std::cout << "2 - best hits calcolati " << std::endl;
    }

    /*
     *
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
        for(auto &map : this->map_best_hits) {
            int id_gene_a = map.first;

            for(auto &submap : map.second) {
                int id_gene_b = submap.first;
                double jaccard = submap.second;

                //it: chiave uguale a gene_b
                auto it = std::find_if(this->map_best_hits.begin(), this->map_best_hits.end(),
                                       [&id_gene_b](const auto& key){ return key.first == id_gene_b;});

                if(it != this->map_best_hits.end()) {

                    //cerca il gene a come chiave della mappa di it (it->second)
                    auto it2 = std::find_if(it->second.begin(), it->second.end(),
                                            [&id_gene_a](const auto& key){ return key.first == id_gene_a;});

                    if(it2 != it->second.end()) {
                        std::unordered_map<int, double> temp;
                        temp.insert(std::make_pair(id_gene_b, jaccard));

                        this->map_bidirectional_best_hits.insert(std::make_pair(id_gene_a, temp));


                        auto it3 = std::find_if(this->vector_tuple_bbh.begin(),
                                                this->vector_tuple_bbh.end(),
                                                [&id_gene_a, &id_gene_b, &jaccard](const std::tuple<int,int,double>& e) {
                                                            return (std::get<0>(e) == id_gene_a &&
                                                                    std::get<1>(e) == id_gene_b &&
                                                                    std::get<2>(e) == jaccard) ||
                                                                    (std::get<0>(e) == id_gene_b &&
                                                                     std::get<1>(e) == id_gene_a &&
                                                                     std::get<2>(e) == jaccard);
                                                });


                        if(it3 == this->vector_tuple_bbh.end())
                            this->vector_tuple_bbh.emplace_back(std::make_tuple(id_gene_a, id_gene_b, jaccard));

                    }
                }
            }
        }

        std::cout << "2 - bidirectional best hits calcolati" << std::endl;
    }

    /*
     * Per ogni occorrenza della mappa, se il gene è diverso, confrontiamo 2 occorrenze alla volta (seconda parte mappa)
     * Cerchiamo ogni occorrenza di questa seconda parte del PRIMO GENE nella seconda parte della mappa del SECONDO GENE
     * Poi lo facciamo al contrario. Se c'è corrispondenza aggiungiamo questo nuovo bbh alla lista
     *
     */
    void calculate_bidirectional_best_hits() {

        /*std::unordered_map<int, int> map_bidirectional_best_hits_internal;

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

                        std::unordered_map<int, double> temp;
                        temp.insert(std::make_pair(d->first, d->second));

                        std::cout << c->first << " " << d->first << std::endl;

                        this->map_bidirectional_best_hits.insert(std::make_pair(c->first, temp));
                    }
                }
            }
        }
        */
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

    std::unordered_map<int, std::unordered_map<int, double>>& get_map_hits() {
        return this->map_hits;
    }

    std::unordered_map<int, std::unordered_map<int, double>>& get_map_best_hits() {
        return this->map_best_hits;
    }

    std::unordered_map<int, std::unordered_map<int, double>>& get_map_bidirectional_best_hits() {
        return this->map_bidirectional_best_hits;
    }

    std::vector<std::tuple<int, int, double>>& get_vector_tuple_bbh() {
        return this->vector_tuple_bbh;
    }

    /*
     * Getter which returns the unsorted set containing the alphabet obtained from all genes in the .faa file
     *
     * @param[out] std::unordered_set<char>&
     */
    std::unordered_set<char>& get_alphabet() {
        return this->alphabet;
    }

private:
    int kmer_size;
    const int flag; //0 amino acids, 1 nucleotides
    std::unordered_set<char> alphabet;
    std::vector<std::vector<int>>* genome_sequenceid;
    std::vector<std::string> sequences;
    std::unordered_map<std::string, std::map<std::bitset<18>, int, bitset_comparer<18>>> sequences_kmers;   //map<sequence - map<kmer, contatore>>
    const std::vector<std::string>* sequences_prefilter;
    std::vector<std::pair<int, int>>* best_hits_prefilter;
    std::unordered_map<int, std::unordered_map<int, double>> map_hits;
    std::unordered_map<int, std::unordered_map<int, double>> map_best_hits;
    std::unordered_map<int, std::unordered_map<int, double>> map_bidirectional_best_hits;
    std::vector<std::tuple<int, int, double>> vector_tuple_bbh;

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

        unsigned long genes_lenght = 0;
        for(auto &i : this->sequences)
            genes_lenght += i.length();

        std::cout << genes_lenght << std::endl;
    }

    void compute_best_hits() {
        for (auto &it: this->map_hits) {
            if(!it.second.empty()) {
                this->map_best_hits.insert(std::make_pair(it.first, std::unordered_map<int, double>()));

                double max_jaccard = 0;

                for(auto &i: it.second) {
                    if(i.second > max_jaccard)
                        max_jaccard = i.second;
                }

                for(auto &i : it.second) {
                    if(i.second == max_jaccard)
                        this->map_best_hits[it.first].insert(std::make_pair(i.first, i.second));
                }
            }
        }
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
