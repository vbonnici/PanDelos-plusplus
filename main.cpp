#include <iostream>
#include <string>
#include "PangeneIData.h"
#include "Omologus.h"
#include "Omologusv2.h"
#include "BestHits.h"
#include "BidirectionalBestHits.h"
#include "PreFilter.h"
#include "Paralog.h"
#include <sys/time.h>
#include "lib/Kvalue.h"
#include "lib/ArgParse.h"

int main(int argc, char* argv[]){
    struct timeval tempo{};
    double t1, t2;

    /*** Argument parsing ***/
    ArgParser parser = ArgParser();
    parser.parse_arguments(argc, argv);

    const char* filename = parser.get_filename();
    int sequences_type = parser.get_sequences_type();
    std::string output = parser.get_output();
    std::string log = parser.get_log();

    std::ofstream log_stream(log, std::ofstream::trunc);
    std::ofstream output_stream(output, std::ofstream::trunc);

    /*** File parsing ***/
        gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    PangeneIData fire = PangeneIData(filename, &log_stream);
    fire.close();

    auto sequences = fire.get_sequences();
    auto genome_sequencesid = fire.get_genome_sequencesid();

    /*** Check kvalue ***/
    Kvalue define_kvalue = Kvalue(&sequences, filename, sequences_type, output, log, &log_stream);
    int kmer_size = define_kvalue.get_kmer_size();
    log_stream << "kmer_size: " << kmer_size << std::endl;

    /*** Prefiltering ***/
    PreFilter filter = PreFilter(sequences, genome_sequencesid, sequences_type, &log_stream);
    filter.init_sequences_kmers();
    filter.calculate_kmer_multiplicity();
    filter.calculate_best_hits();
    auto prefilter_best_hits = filter.get_best_hits();

    /*** Omologus ***/
    Omologusv2 omologus = Omologusv2(sequences, prefilter_best_hits, sequences_type, kmer_size, &log_stream);
    omologus.init_sequences_kmers();
    omologus.calculate_kmer_multiplicity();
    omologus.calculate_best_hits();
    auto map_hits = omologus.get_map_hits();

    /*** BestHits ***/
    BestHits bh = BestHits(map_hits, &log_stream);
    bh.compute_best_hits();
    //Helper::nested_unordered_map_print<int, int, double>(output_stream, map_hits, " ");
    auto map_best_hits = bh.get_map_best_hits();

    /*** BidirectionalBestHits ***/
    BidirectionalBestHits bbh = BidirectionalBestHits(map_hits, &log_stream); //map_best_hits
    bbh.calculate_bbh();
    auto vector_tuple_bbh = bbh.get_vector_tuple_bbh();

    /*** Paralog ***/
    Paralog paralog = Paralog(sequences, genome_sequencesid, sequences_type, kmer_size, vector_tuple_bbh, &log_stream);
    paralog.calculate_paralog();
    auto paralog_best_hits = paralog.get_paralog_best_hits();


    /*** Output ***/
    output_stream << "Ortologhi " << std::endl;

    for(auto &i : vector_tuple_bbh)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");

    output_stream << "Paraloghi: " << std::endl;

    for(auto &i : paralog_best_hits)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");


    output_stream.close();

        gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    log_stream << "Tempo computazione: " << t2-t1 << std::endl;

    return 0;
};