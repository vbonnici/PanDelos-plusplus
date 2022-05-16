#include "PangeneIData.h"
#include "PreFiltering.h"
#include <sys/time.h>

int main(int argc, char* argv[]){

    struct timeval tempo{};
    double t1, t2;

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(argv[1]);

    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    fire.close();
    fire.compute_alphabet();

    std::vector<int> sequences_genome = fire.get_sequences_genome();
    std::vector<std::string> genomes_names = fire.get_genomes_names();
    std::vector<std::string> sequences_name = fire.get_sequences_name();
    std::vector<std::string> sequences_description = fire.get_sequences_description();
    std::vector<std::string> sequences = fire.get_sequences();
    std::unordered_set<char> alphabet = fire.get_alphabet();

    PreFiltering filter = PreFiltering(6, sequences);

    filter.populate_map_sequences(); //tempo di esecuzione troppo lungo

    filter.calculate_kmer_frequency();

    auto result = filter.get_map_sequences();


    std::cout << "Tempo di lettura del file (secondi): " << t2-t1 << std::endl;

    return 0;
};