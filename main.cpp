#include "PangeneIData.h"
#include <sys/time.h>
using namespace std;

int main(int argc, char* argv[]){

    struct timeval tempo{};
    double t1, t2;

    gettimeofday(&tempo,nullptr); t1 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    PangeneIData fire = PangeneIData(argv[1]);

    gettimeofday(&tempo,nullptr); t2 = tempo.tv_sec+(tempo.tv_usec/1000000.0);

    fire.close();
    fire.compute_alphabet();

    //fire.print_genomes_names();
    //fire.print_sequences_name();
    //fire.print_sequences();
    //fire.print_sequences_description();
    //fire.print_alphabet();

    //fire.print_sequences_genome();

    cout << "Tempo di lettura del file (secondi): " << t2-t1 << endl;

    return 0;
};