#include <iostream>
#include <string>
#include "../../include/PangeneIData.h"
#include "../../lib/PanDelos/Homologues.h"
#include "../../lib/PanDelos/BestHits.h"
#include "../../lib/PanDelos/BidirectionalBestHits.h"
#include "../../lib/PanDelos/PreFilter.h"
#include "../../lib/PanDelos/Paralogues.h"
#include <sys/time.h>
#include "../../include/Kvalue.h"
#include "../../lib/Argparser/ArgParse.h"

int main(int argc, char* argv[]){
    struct timeval time{};
    double total_time_start, total_time_end;
    double prefiltering_start, prefiltering_end;
    double homologues_start, homologues_end;
    double bh_start, bh_end;
    double bbh_start, bbh_end;
    double paralogues_start, paralogues_end;

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
        gettimeofday(&time,nullptr); total_time_start = time.tv_sec+(time.tv_usec/1000000.0);
    PangeneIData file_reader = PangeneIData(filename, &log_stream);
    file_reader.close();

    auto sequences = file_reader.get_sequences();
    auto genome_sequencesid = file_reader.get_genome_sequencesid();
    auto genes_id_interval = file_reader.get_genes_id_interval();

    /*** Check kvalue ***/
    Kvalue define_kvalue = Kvalue(&sequences, filename, sequences_type, output, log, &log_stream);
    int kmer_size = define_kvalue.get_kmer_size();

    /*** Prefiltering ***/
    std::vector<std::pair<int, int>> pre_filter_candidate_sequences;
    try {
            gettimeofday(&time,nullptr); prefiltering_start = time.tv_sec+(time.tv_usec/1000000.0);
        PreFilter pre_filter = PreFilter(sequences, genome_sequencesid, sequences_type, &log_stream);
        pre_filter.find_candidate_sequences();
        pre_filter_candidate_sequences = pre_filter.get_candidate_sequences();
            gettimeofday(&time,nullptr); prefiltering_end = time.tv_sec+(time.tv_usec/1000000.0);
        log_stream << "Prefilter phase completed in " << prefiltering_end-prefiltering_start << " seconds" << std::endl;

    } catch(std::exception const& e) {
        log_stream << "Exception: " << e.what() << "\n";
        exit(11);
    }


    /*** Homologues ***/
    std::unordered_map<int, std::unordered_map<int, double>> homologues_candidate_sequences;
    try {
            gettimeofday(&time,nullptr); homologues_start = time.tv_sec+(time.tv_usec/1000000.0);
        Homologues homologues = Homologues(sequences, pre_filter_candidate_sequences, sequences_type, kmer_size, &log_stream);
        homologues.find_candidate_sequences();
        homologues_candidate_sequences = homologues.get_candidate_sequences();
            gettimeofday(&time,nullptr); homologues_end = time.tv_sec+(time.tv_usec/1000000.0);
        log_stream << "Homologues phase completed in " << homologues_end-homologues_start << " seconds" << std::endl;

    } catch(std::exception const& e) {
        log_stream << "Exception: " << e.what() << "\n";
        exit(11);
    }

    /*** BestHits ***/
    std::unordered_map<int, std::unordered_map<int, double>> best_hits;
    try {
            gettimeofday(&time,nullptr); bh_start = time.tv_sec+(time.tv_usec/1000000.0);
        BestHits bh = BestHits(homologues_candidate_sequences, genes_id_interval, &log_stream);
        bh.find_best_hits();
        best_hits = bh.get_best_hits();
            gettimeofday(&time,nullptr); bh_end = time.tv_sec+(time.tv_usec/1000000.0);
        log_stream << "Best Hits phase completed in " << bh_end-bh_start << " seconds" << std::endl;

    } catch(std::exception const& e) {
        log_stream << "Exception: " << e.what() << "\n";
        exit(11);
    }


    /*** BidirectionalBestHits ***/
    std::vector<std::tuple<int, int, double>> bidirectional_best_hits;
    try {
            gettimeofday(&time,nullptr); bbh_start = time.tv_sec+(time.tv_usec/1000000.0);
        BidirectionalBestHits bbh = BidirectionalBestHits(best_hits, &log_stream);
        bbh.find_bidirectional_best_hits();
        bidirectional_best_hits = bbh.get_bidirectional_best_hits();
            gettimeofday(&time,nullptr); bbh_end = time.tv_sec+(time.tv_usec/1000000.0);
        log_stream << "Bidirectional Best Hits phase completed in " << bbh_end-bbh_start << " seconds" << std::endl;

    } catch(std::exception const& e) {
        log_stream << "Exception: " << e.what() << "\n";
        exit(11);
    }


    /*** Paralogues ***/
    std::vector<std::tuple<int, int, double>> paralogues_best_hits;
    try {
            gettimeofday(&time,nullptr); paralogues_start = time.tv_sec+(time.tv_usec/1000000.0);
        Paralogues paralogues = Paralogues(sequences, genome_sequencesid, genes_id_interval, sequences_type, kmer_size, bidirectional_best_hits, &log_stream);
        paralogues.find_paralogues();
        paralogues_best_hits = paralogues.get_paralogues_best_hits();
            gettimeofday(&time,nullptr); paralogues_end = time.tv_sec+(time.tv_usec/1000000.0);
        log_stream << "Paralogues phase completed in " << paralogues_end-paralogues_start << " seconds" << std::endl;

    } catch(std::exception const& e) {
        log_stream << "Exception: " << e.what() << "\n";
        exit(11);
    }


    /*** Output ***/
    //output_stream << "Orthologues " << std::endl;

    for(auto &i : bidirectional_best_hits)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");

    //output_stream << "Paralogues: " << std::endl;

    for(auto &i : paralogues_best_hits)
        Helper::simple_triple_print<int, int, double>(output_stream, i, "\t");


    output_stream.close();

        gettimeofday(&time,nullptr); total_time_end = time.tv_sec+(time.tv_usec/1000000.0);

    log_stream << "Computation time net of preliminary operations: " << total_time_end-total_time_start << " seconds" << std::endl;

    return 0;
};