#ifndef PANDELOS_PLUSPLUS_PANGENEIDATA_H
#define PANDELOS_PLUSPLUS_PANGENEIDATA_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

class PangeneIData {
public:
    explicit PangeneIData(const char* filename) {
        this->open_file(filename);

        bool nameline = true;
        std::string genome_name, sequence_name, sequence_description;
        std::map<std::string, int> genomeID;
        int genomeid;

        while (this->is_valid()) {
            if (nameline) {
                genome_name = this->next_string();
                sequence_name = this->next_string();
                sequence_description = this->next_string();
            } else {
                this->sequences.emplace_back(sequence_name + "," + this->next_string());
                this->sequences_name.push_back(sequence_name);
                auto genomeid_iterator = genomeID.find(genome_name);

                if (genomeid_iterator == genomeID.end()) {
                    genomeid = genomeID.size();
                    genomeID.insert(std::pair<std::string, int>(genome_name, genomeid));
                }
                else
                    genomeid = genomeid_iterator->second;

                this->sequences_genome.push_back(genomeid);
                this->sequences_description.push_back(sequence_description);
            }

            nameline = !nameline;
        }

        this->genomes_names.reserve(genomeID.size());
        for (auto &it: genomeID)
            this->genomes_names.push_back(it.first);
    }

    void compute_alphabet() {
        for (auto &i: this->sequences)
            for(char &a : i)
                this->alphabet.insert(a);
    }


    void print_genomes_names() {
        print_container<std::string, std::vector<std::string>::iterator>
            (std::cout, this->genomes_names.begin(), this->genomes_names.end(), "\n");
    }

    void print_sequences_name() {
        print_container<std::string, std::vector<std::string>::iterator>
                (std::cout, this->sequences_name.begin(), this->sequences_name.end(), "\n");
    }

    void print_sequences_description() {
        print_container<std::string, std::vector<std::string>::iterator>
                (std::cout, this->sequences_description.begin(), this->sequences_description.end(), "\n");
    }

    void print_sequences() {
        print_container<std::string, std::vector<std::string>::iterator>
                (std::cout, this->sequences.begin(), this->sequences.end(), "\n");
    }

    void print_alphabet() {
        print_container<char, std::unordered_set<char>::iterator>
                (std::cout, this->alphabet.begin(), this->alphabet.end(), "\n");
    }

    void print_sequences_genome() {
        print_container<int, std::vector<int>::iterator>
                (std::cout, this->sequences_genome.begin(), this->sequences_genome.end(), "\n");
    }

    std::vector<int>& get_sequences_genome() {
        return this->sequences_genome;
    }

    std::vector<std::string>& get_genomes_names() {
        return this->genomes_names;
    }

    std::vector<std::string>& get_sequences_name() {
        return this->sequences_name;
    }

    std::vector<std::string>& get_sequences_description() {
        return this->sequences_description;
    }

    std::vector<std::string>& get_sequences() {
        return this->sequences;
    }

    std::unordered_set<char>& get_alphabet() {
        return this->alphabet;
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

    std::vector<std::string> sequences;
    std::vector<std::string> sequences_name; //TODO: forse non serve
    std::vector<std::string> sequences_description;
    std::vector<int> sequences_genome;
    std::vector<std::string> genomes_names; //TODO: forse non serve
    std::unordered_set<char> alphabet;

    //TODO: aggiungere exception safety apertura file
    void open_file(const char* filename) {
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

    template<typename T, typename InputIterator>
    void print_container(std::ostream& ostr, InputIterator itbegin, InputIterator itend, const std::string& delimiter) {
        std::copy(itbegin, itend, std::ostream_iterator<T>(ostr, delimiter.c_str()));
    }
};

#endif //PANDELOS_PLUSPLUS_PANGENEIDATA_H