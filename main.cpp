#include "PangeneIData.h"
#include "BidirectionalBestHits.h"
#include "BidirectionalBestHitsV2.h"
#include "PreFilter.h"
#include <sys/time.h>

int main(int argc, char* argv[]){

    struct timeval tempo{};
    double t1, t2;

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(argv[1]);

    fire.close();

    auto sequences = fire.get_sequences();
    auto map_sequences_attributes = fire.get_map_sequences_attributes();
    auto sequences_name = fire.get_sequences_name();
    auto genome_sequencesid = fire.get_genome_sequencesid();

    PreFilter filter = PreFilter(sequences, genome_sequencesid, 1);

    filter.init_sequences_kmers();

    filter.calculate_kmer_multiplicity();

    filter.calculate_best_hits_v2();

    auto best_hits = filter.get_best_hits();

    BidirectionalBestHitsV2 bbh = BidirectionalBestHitsV2(sequences, best_hits, 1);

    bbh.init_sequences_kmers();

    bbh.calculate_kmer_multiplicity();

    auto map_sequences_kmers = bbh.get_sequences_kmers();

    bbh.calculate_best_hits();

    auto map_best_hits_2 = bbh.get_map_best_hits();
/*
    bbh.calculate_bidirectional_best_hits();

    auto map_bidirectional_best_hits = bbh.get_map_bidirectional_best_hits();

    for(auto &i : map_bidirectional_best_hits) {
        auto second_gene = i.second;
        auto second_gene_value = second_gene.begin();

        std::cout << i.first << " " << second_gene_value->first << " " << second_gene_value->second << std::endl;
    }
*/
/*    BidirectionalBestHits bbh = BidirectionalBestHits(sequences, map_best_hits, map_sequences_attributes, 1);

    bbh.init_map_sequences_kmers();

    bbh.calculate_kmer_multiplicity();

    auto map_sequences_kmers = bbh.get_map_sequences_kmers();

    bbh.calculate_best_hits();

    auto map_best_hits_2 = bbh.get_map_best_hits();

    bbh.calculate_bidirectional_best_hits();

    auto map_bidirectional_best_hits = bbh.get_map_bidirectional_best_hits();

    for(auto &i : map_bidirectional_best_hits) {
        auto second_gene = i.second;
        auto second_gene_value = second_gene.begin();

        std::cout << i.first << " " << second_gene_value->first << " " << second_gene_value->second << std::endl;
    }
*/
    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    std::cout << "Tempo computazione: " << t2-t1 << std::endl;

    return 0;
};