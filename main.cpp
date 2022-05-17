#include "PangeneIData.h"
#include "PreFiltering.h"
#include <sys/time.h>

int main(int argc, char* argv[]){

    struct timeval tempo{};
    double t1, t2;

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(argv[1]);

    //gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    fire.close();
    fire.compute_alphabet();

    auto sequences_genome = fire.get_sequences_genome();
    auto genomes_names = fire.get_genomes_names();
    auto sequences_name = fire.get_sequences_name();
    auto sequences_description = fire.get_sequences_description();
    auto sequences = fire.get_sequences();
    auto alphabet = fire.get_alphabet();

    PreFiltering filter = PreFiltering(6, sequences);

    filter.populate_map_sequences();

    filter.calculate_kmer_frequency();

    auto map_sequences = filter.get_map_sequences();

    filter.calculate_bh();

    auto map_genes_jaccard = filter.get_map_genes_jaccard();

    filter.calculate_bbh();

    auto genes_bbh = filter.get_genes_bbh();

    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    std::cout << "Tempo computazione: " << t2-t1 << std::endl;

    return 0;
};