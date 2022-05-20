#include "PangeneIData.h"
#include "PreFiltering.h"
#include "BestHits.h"
#include <sys/time.h>

int main(int argc, char* argv[]){

    struct timeval tempo{};
    double t1, t2;

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(argv[1]);

    //gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    fire.close();
    //fire.compute_alphabet();

    //auto sequences_genome = fire.get_sequences_genome();
    //auto genomes_names = fire.get_genomes_names();
    //auto sequences_name = fire.get_sequences_name();
    //auto sequences_description = fire.get_sequences_description();
    auto sequences = fire.get_sequences();
    //auto alphabet = fire.get_alphabet();

    BestHits bh = BestHits(sequences, 0);

    bh.init_map_sequences_kmers();

    auto map_sequences_kmers_before = bh.get_map_sequences_kmers();

    bh.calculate_kmer_multiplicity_nucleotides();

    auto map_sequences_kmers_after = bh.get_map_sequences_kmers();

    //std::cout << map_sequences_kmers_after << std::endl;

    //std::cout << map_sequences_kmers_after[0]["00_genome:NZ_CP007390.1:CF57_RS00010:11"][9] << std::endl;

    /*std::unordered_map<std::string, int*> map_sequences_kmers;

    int* array = new int[4095];
    std::fill_n(array, 4095, 0);

    map_sequences_kmers.insert(std::make_pair("ciao", array));

    delete[] array;

    map_sequences_kmers["ciao"][9] = 9;
     */

    /*std::unordered_map<std::string, int*> test;

    int array[4095];
    std::fill_n(array, 4095, 1);

    test.insert(std::make_pair("ciao", array));

    test["ciao"][2] = 9;

    std::cout << test["ciao"][4000] << std::endl;
    */
    /*PreFiltering filter = PreFiltering(6, sequences);

    filter.init_map_sequences_kmers();

    filter.calculate_kmer_frequency();

    auto map_sequences_kmers = filter.get_map_sequences_kmers();

    filter.calculate_best_hits();

    auto map_best_hits = filter.get_map_best_hits();

    filter.calculate_bidirectional_best_hits();

    auto map_bidirectional_best_hits = filter.get_map_bidirectional_best_hits();

     */
    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    std::cout << "Tempo computazione: " << t2-t1 << std::endl;

    return 0;
};