#include <iostream>
#include <string>
#include "PangeneIData.h"
#include "BidirectionalBestHits.h"
#include "PreFilter.h"
#include <sys/time.h>
#include "lib/Kvalue.h"
#include "lib/helper.h"


int main(int argc, char* argv[]){
    Helper helper = Helper();

    helper.parse_arguments(argc, argv);

    int kmer_size;

    struct timeval tempo{};
    double t1, t2;

    const char* filename = helper.get_filename();
    int sequences_type = helper.get_sequences_type();
    std::string output = helper.get_output();
    std::string log = helper.get_log();

    std::ofstream log_stream(log, std::ofstream::trunc);

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(filename);

    fire.close();

    auto sequences = fire.get_sequences();
    auto genome_sequencesid = fire.get_genome_sequencesid();

    Kvalue define_kvalue = Kvalue(&sequences, filename, sequences_type, output, log);
    kmer_size = define_kvalue.get_kmer_size();

    std::cout << kmer_size << std::endl;

    PreFilter filter = PreFilter(sequences, genome_sequencesid, sequences_type);

    filter.init_sequences_kmers();

    filter.calculate_kmer_multiplicity();

    filter.calculate_best_hits();

    auto prefilter_best_hits = filter.get_best_hits();

    BidirectionalBestHits bbh = BidirectionalBestHits(sequences, prefilter_best_hits, genome_sequencesid, sequences_type);

    bbh.init_sequences_kmers();

    bbh.calculate_kmer_multiplicity();

    bbh.calculate_best_hits();

    auto map_best_hits = bbh.get_map_best_hits();

    bbh.calculate_bbh();

    auto vector_tuple_bbh = bbh.get_vector_tuple_bbh();

    std::ofstream output_stream(output, std::ofstream::trunc);

    for(auto &i : vector_tuple_bbh) {
        output_stream << std::get<0>(i) << " " << std::get<1>(i) << " " << std::get<2>(i) << std::endl;
    }

    output_stream.close();

    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    log_stream << "Tempo computazione: " << t2-t1 << std::endl;

    return 0;
};