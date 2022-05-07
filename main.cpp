#include "PangeneIData.h"

int main(int argc, char* argv[]){

    std::cout<<argv[1]<<"\n";

    PangeneIData fire = PangeneIData(argv[1]);

    /* DEBUG VECTOR
     * fire.print_sequences();
    fire.print_sequenceName();
    fire.print_sequenceDescription();
    fire.get_genomeNames();
     */

    fire.close();


    return 0;
};