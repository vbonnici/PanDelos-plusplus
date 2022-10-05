#ifndef PANDELOS_PLUSPLUS_HOMOLOGUES
#define PANDELOS_PLUSPLUS_HOMOLOGUES

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <algorithm>
#include "../../conf/conf.h"
#include <omp.h>
#include "../../include/Helper.h"

/*
 * Class that calculates the similarity between gene sequences, being able to use a medium-high kmer size value.
 */
class Homologues {
public:

    /*
     * The constructor constructs the adequate data structures in order to proceed with the analysis and takes care of
     * calling the methods that enumerate the k-mer and calculate the multiplicities
     *
     * @param[in] const std::vector<std::string>&
     * @param[in] const std::vector<std::pair<int, int>>&
     * @param[in] const int
     * @param[in] const int
     * @param[in] std::ofstream*
     *
     * @throws std::runtime_error if the log_stream's badbit error state flag is set
     * @throws std::runtime_error if the input sequences vector is empty
     * @throws std::runtime_error if the genome sequencesid vector is empty
     */
    explicit Homologues(const std::vector<std::string>& input_sequences,
                        const std::vector<std::pair<int, int>>& gene_pair_input,
                        const int sequences_type,
                        const int kmer_size,
                        std::ofstream* log_stream) : sequences_type(sequences_type), kmer_size(kmer_size){

        if(log_stream->bad())
            throw std::runtime_error("the log_stream's badbit error state flag is set");

        if(input_sequences.empty())
            throw std::runtime_error("input sequences vector is empty");

        if(gene_pair_input.empty())
            throw std::runtime_error("gene_pair_input vector is empty");

        this->log_stream = log_stream;
        this->input_sequences = &input_sequences;
        this->gene_pair_input = &gene_pair_input;

        this->collect_sequences_id_from_gene_pair(&gene_pair_input);
        this->init_sequences_kmers();
        this->calculate_kmer_multiplicity();
    }

    /*
     * Main method of the class that deals with finding candidate sequences to be homologous.
     *
     * The method works on a data structure containing identifiers of gene pairs to be analyzed received as input at
     * the time of construction.
     *
     * For each gene of each pair to be analyzed, a data structure is created in which the identifiers of the genes
     * most similar to it will converge.
     *
     * The homologues phase consists, in parallel, of comparing the multiplicity of k-mer for pairs of genomes
     * in order to calculate the Jaccard similarity index, which establishes the similarity between two genes.
     *
     * This comparison is obtained by combining in a single data structure the k-mer found by the 2 genes as well as
     * relative multiplicities
     *
     * @param[in] double
     */
    void find_candidate_sequences(const double jaccard_threshold = 0.0) {
        int id_gene_a;
        int id_gene_b;
        unsigned int value_a;
        unsigned int value_b;
        unsigned int counter_min;
        unsigned int counter_max;
        double jaccard_similarity;

        for (auto &i: *this->gene_pair_input) {
            this->candidate_sequences.insert(std::make_pair(i.first, std::unordered_map<int, double>()));
            this->candidate_sequences.insert(std::make_pair(i.second, std::unordered_map<int, double>()));
        }

        omp_set_num_threads(omp_get_num_procs());

        #pragma omp parallel for private(jaccard_similarity, counter_min, counter_max, value_a, value_b, id_gene_a, id_gene_b)
        for (auto &i: *this->gene_pair_input) {
            id_gene_a = i.first;
            id_gene_b = i.second;

            counter_min = 0;
            counter_max = 0;

            auto kmer_key_a = this->sequences_kmers.find(id_gene_a);
            auto kmer_key_b = this->sequences_kmers.find(id_gene_b);

            std::map<std::string, int> kmerstring_map_a;

            auto kmerbit_map_a = kmer_key_a->second;

            for(auto &kmerbit : kmerbit_map_a)
                kmerstring_map_a.insert(std::make_pair(kmerbit.first.to_string(), kmerbit.second));

            std::map<std::string, int> kmerstring_map_b;

            auto kmerbit_map_b = kmer_key_b->second;

            for(auto &kmerbit : kmerbit_map_b)
                kmerstring_map_b.insert(std::make_pair(kmerbit.first.to_string(), kmerbit.second));

            auto kmerstring_map_a_b = Helper::UnionMaps<std::string, int, int>(kmerstring_map_a, kmerstring_map_b);

            for(auto &kmer : kmerstring_map_a_b) {

                auto kmercounter_a_b = kmer.second;

                value_a = kmercounter_a_b.first;
                value_b = kmercounter_a_b.second;

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

            if (counter_max > 0 && jaccard_similarity >= jaccard_threshold) {
                auto result_a = this->candidate_sequences.find(id_gene_a);

                #pragma omp critical
                result_a->second.insert(std::make_pair(id_gene_b, jaccard_similarity));


                auto result_b = this->candidate_sequences.find(id_gene_b);

                #pragma omp critical
                result_b->second.insert(std::make_pair(id_gene_a, jaccard_similarity));
            }
        }
    }

    /*
     * Getter which returns the vector of k-mer extracted for each sequence
     *
     * @param[out] const std::vector<int>&
     * @throws std::runtime_error if sequences id vector is empty
     */
    const std::vector<int>& get_sequences_id() {
        if(this->sequences_id.empty())
            throw std::runtime_error("sequences id vector is empty");

        return this->sequences_id;
    }

    /*
     * Getter which returns the set of k-mer extracted for each sequence
     *
     * @param[out] const std::unordered_map<int, std::map<std::bitset<kvalue>, int, Helper::bitset_comparer<kvalue>>>&
     * @throws std::runtime_error if sequences kmers vector is empty
     */
    const std::unordered_map<int, std::map<std::bitset<kvalue>, int, Helper::bitset_comparer<kvalue>>>& get_sequences_kmers() {
        if(this->sequences_kmers.empty())
            throw std::runtime_error("sequences kmers vector is empty");

        return this->sequences_kmers;
    }

    /*
     * Getter that returns the set of candidate sequences that have passed the homologues phase
     *
     * @param[out] const std::vector<std::pair<int, int>>&
     */
    const std::unordered_map<int, std::unordered_map<int, double>>& get_candidate_sequences() {
        return this->candidate_sequences;
    }


private:
    std::ofstream* log_stream;
    int kmer_size;
    const int sequences_type;
    std::vector<int> sequences_id;
    std::unordered_map<int, std::map<std::bitset<kvalue>, int, Helper::bitset_comparer<kvalue>>> sequences_kmers;   //map<sequence_id - map<kmer, contatore>>
    const std::vector<std::string>* input_sequences;
    const std::vector<std::pair<int, int>>* gene_pair_input;
    std::unordered_map<int, std::unordered_map<int, double>> candidate_sequences;

    /*
     * This method sequentially inserts pairs of genes present in a vector of pair into a single vector
     *
     * @param[in] const std::vector<std::pair<int, int>>*
     */
    void collect_sequences_id_from_gene_pair(const std::vector<std::pair<int, int>>* gene_pair_in) {

        std::set<int> temp_sequences;

        for(auto &i : *gene_pair_in) {
            temp_sequences.insert(i.first);
            temp_sequences.insert(i.second);
        }

        this->sequences_id.assign(temp_sequences.begin(), temp_sequences.end());
    }

    /*
     * Initializes a map for each sequence with the first k-mer and the counter at 0
     */
    void init_sequences_kmers() {

        for(auto &sequence_id : this->sequences_id) {

            std::map<std::bitset<kvalue>, int, Helper::bitset_comparer<kvalue>> temp;
            std::bitset<kvalue> bitset_temp;
            bitset_temp.set(false);
            temp.insert(std::make_pair(bitset_temp, 0));

            this->sequences_kmers.insert(std::make_pair(sequence_id, temp));
        }
    }

    /*
     *
     * The multiplicity of k-mer is calculated with a sliding window algorithm which takes care of traversing
     * the entire sequence for all possible substrings of calculated kmer size length.
     * The extracted k-mer are transformed into bits and each corresponds to a unique index of the previously created array.
     * Then the enumeration takes place by incrementing the counters for certain bits from time to time
     */
    void calculate_kmer_multiplicity() {
        std::string kmer;
        for(auto &i : this->sequences_kmers) {
            std::string sequence = this->input_sequences->operator[](i.first);

            for(int window = 0; window < sequence.length(); ++window) {
                if(this->sequences_type == 1)
                    kmer = sequence.substr(window, this->kmer_size);
                else {
                    std::string aminoacid = sequence.substr(window, this->kmer_size);
                    kmer = Helper::aminoacid_to_nucleotides(aminoacid);
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
    }

    /*
     * Method that checks if a given kmer extracted is valid
     *
     * @param[in] const std::string
     * @param[out] bool
     */
    [[nodiscard]] bool kmer_is_valid(const std::string &str) const {
        if(this->sequences_type == 0)
            return str.length() == this->kmer_size*3 && str.find_first_not_of("ACGT") == std::string::npos;
        else
            return str.length() == this->kmer_size && str.find_first_not_of("ACGT") == std::string::npos;
    }

    /*
     * Method that converts each k-mer into a sequence of bits using bitwise operations
     *
     * The sequence is stored in the bitset data structure
     *
     * @param[in] const std::string&
     * @param[out] std::bitset<kvalue>
     */
    static std::bitset<kvalue> kmer_to_bit(const std::string& kmer) {
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
};
#endif //PANDELOS_PLUSPLUS_HOMOLOGUES