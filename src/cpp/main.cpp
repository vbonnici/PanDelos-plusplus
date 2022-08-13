#include <iostream>
#include <string>
#include "../../include/PangeneIData.h"
#include "../../lib/PanDelos/Omologus.h"
#include "../../lib/PanDelos/BestHits.h"
#include "../../lib/PanDelos/BidirectionalBestHits.h"
#include "../../lib/PanDelos/PreFilter.h"
#include "../../lib/PanDelos/Paralog.h"
#include <sys/time.h>
#include "../../include/Kvalue.h"
#include "../../lib/Argparser/ArgParse.h"

int main(int argc, char* argv[]){
    struct timeval tempo{};
    double net_start, net_end;
    double prefiltering_start, prefiltering_end;
    double omologus_start, omologus_end;
    double bh_start, bh_end;
    double bbh_start, bbh_end;
    double paralog_start, paralog_end;

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
        gettimeofday(&tempo,nullptr); net_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    PangeneIData fire = PangeneIData(filename, &log_stream);
    fire.close();

    auto sequences = fire.get_sequences();
    auto genome_sequencesid = fire.get_genome_sequencesid();
    auto genes_id_interval = fire.get_genes_id_interval();

    /*** Check kvalue ***/
    Kvalue define_kvalue = Kvalue(&sequences, filename, sequences_type, output, log, &log_stream);
    int kmer_size = define_kvalue.get_kmer_size();

    /*** Prefiltering ***/
        gettimeofday(&tempo,nullptr); prefiltering_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    PreFilter filter = PreFilter(sequences, genome_sequencesid, sequences_type, &log_stream);
    filter.init_sequences_kmers();
    filter.calculate_kmer_multiplicity();
    filter.calculate_best_hits();
    auto prefilter_best_hits = filter.get_best_hits();
        gettimeofday(&tempo,nullptr); prefiltering_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    log_stream << "Prefilter phase completed in " << prefiltering_end-prefiltering_start << " seconds" << std::endl;

    /*** Omologus ***/
        gettimeofday(&tempo,nullptr); omologus_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    Omologus omologus = Omologus(sequences, prefilter_best_hits, sequences_type, kmer_size, &log_stream);
    omologus.init_sequences_kmers();
    omologus.calculate_kmer_multiplicity();
    omologus.calculate_best_hits();
    auto map_hits = omologus.get_map_hits();
        gettimeofday(&tempo,nullptr); omologus_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    log_stream << "Omologus phase completed in " << omologus_end-omologus_start << " seconds" << std::endl;

    /*** BestHits ***/
        gettimeofday(&tempo,nullptr); bh_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    BestHits bh = BestHits(map_hits, genes_id_interval, &log_stream);
    bh.compute_best_hits();
    auto map_best_hits = bh.get_map_best_hits();
        gettimeofday(&tempo,nullptr); bh_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    log_stream << "Best Hits calculated in " << bh_end-bh_start << " seconds" << std::endl;


    /*** BidirectionalBestHits ***/
        gettimeofday(&tempo,nullptr); bbh_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    BidirectionalBestHits bbh = BidirectionalBestHits(map_best_hits, &log_stream);
    bbh.calculate_bbh();
    auto vector_tuple_bbh = bbh.get_vector_tuple_bbh();
        gettimeofday(&tempo,nullptr); bbh_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    log_stream << "Bidirectional Best Hits calculated in " << bbh_end-bbh_start << " seconds" << std::endl;

    /*** Paralog ***/
        gettimeofday(&tempo,nullptr); paralog_start = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    Paralog paralog = Paralog(sequences, genome_sequencesid, genes_id_interval, sequences_type, kmer_size, vector_tuple_bbh, &log_stream);
    paralog.calculate_paralog();
    auto paralog_best_hits = paralog.get_paralog_best_hits();
        gettimeofday(&tempo,nullptr); paralog_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);
    log_stream << "Paralogs calculated in " << paralog_end-paralog_start << " seconds" << std::endl;


    /*** Output ***/
    //output_stream << "Orthologues " << std::endl;

    for(auto &i : vector_tuple_bbh)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");

    //output_stream << "Paralogues: " << std::endl;

    for(auto &i : paralog_best_hits)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");


    output_stream.close();

        gettimeofday(&tempo,nullptr); net_end = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    log_stream << "Computation time net of preliminary operations: " << net_end-net_start << " seconds" << std::endl;

    return 0;
};