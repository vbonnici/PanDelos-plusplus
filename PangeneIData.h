#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <map>
using namespace std;

class PangeneIData {
public:
    explicit PangeneIData(const char *filename) {
        this->init(filename);

        string line;
        bool nameLine = true;
        string genomeName, seqName, product, seq;
        vector<string> cc;
        map<string, int> genomeID;
        int genomeid;

        while (this->is_valid()) {
            if (nameLine) {
                genomeName = this->next_string();
                seqName = this->next_string();
                product = this->next_string();
            } else {
                this->sequences.emplace_back(this->next_string());
                this->sequenceName.push_back(seqName);
                auto genomeid_iterator = genomeID.find(genomeName);

                if (genomeid_iterator == genomeID.end()) {
                    genomeid = genomeID.size();
                    genomeID.insert(pair<string, int>(genomeName, genomeid));
                }
                else
                    genomeid = genomeid_iterator->second;

                this->sequenceGenome.push_back(genomeid);
                this->sequenceDescription.push_back(product);
            }

            nameLine = !nameLine;
        }

        this->genomeNames.reserve(genomeID.size());
        for (auto &it: genomeID)
            this->genomeNames.push_back(it.first);
    }

    //TODO: funzione templatica print vector
    void print_sequences() {
        for (auto & i: this->sequences)
            std::cout << i << endl;
    }

    void print_sequenceName() {
        for (auto & i: this->sequenceName)
            std::cout << i << endl;
    }

    void print_sequenceDescription() {
        for (auto & i: this->sequenceDescription)
            std::cout << i << endl;
    }

    void get_genomeNames() {
        for (auto & i: this->sequenceName)
            std::cout << i << endl;
    }

   void close() {
        fclose(this->pFile);
        free(this->buffer);
    }

private:
    FILE* pFile{};
    long lSize{};
    char* buffer{};
    long pi{};

    vector<string> sequences;
    vector<string> sequenceName;
    vector<string> sequenceDescription;
    vector<int> sequenceGenome;
    vector<string> genomeNames;

    //TODO: aggiungere exception safety apertura file
    void init(const char *filename) {
        size_t result;

        pFile = fopen(filename, "rb");
        if (pFile == NULL) {
            fputs("File error\n", stderr);
            exit(1);
        }

        // obtain file size:
        fseek(pFile, 0, SEEK_END);
        lSize = ftell(pFile);
        rewind(pFile);

        // allocate memory to contain the whole file:
        buffer = (char*) malloc(sizeof(char) * lSize);
        if (buffer == NULL) {
            fputs("Memory error\n", stderr);
            exit(2);
        }

        // copy the file into the buffer:
        result = fread(buffer, 1, lSize, pFile);
        if (result != lSize) {
            fputs("Reading error\n", stderr);
            exit(3);
        }

        pi = 0;

        //std::cout<<buffer<<"\n";
    }

    [[nodiscard]] bool is_valid() const {
        return this->pi < this->lSize;
    }

    char* next_string() {

        while ((this->is_valid()) && ((buffer[pi] == '\n') || (buffer[pi] == '\t') || (buffer[pi] == '\r'))) {
            pi++;
        }

        long ci = pi;

        while ((this->is_valid()) && (buffer[pi] != '\n') && (buffer[pi] != '\t') && (buffer[pi] != '\r')) {
            pi++;
        }

        if (this->is_valid()) {
            buffer[pi] = '\0';
        }
        pi++;

        return buffer + ci;
    }
};