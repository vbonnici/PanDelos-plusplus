#ifndef PANDELOS_PLUSPLUS_KVALUE_H
#define PANDELOS_PLUSPLUS_KVALUE_H

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include "../include/kvalue.h"
#include <stdlib.h>

class Kvalue {
public:

    /*
     * For each gene sequence, it inserts each character of which it is composed into an unordered set
     * (if not yet present)
     */
    explicit Kvalue(const std::vector<std::string>* sequences_input, const char* filename, int sequences_type, std::string output, std::string logfile) {
        for (auto &i: *sequences_input) {
            for(char a : i)
                this->alphabet.insert(a);
        }

        std::cout << "alfabeto calcolato" << std::endl;

        this->genes_lenght = 0;
        for(auto &i : *sequences_input)
            this->genes_lenght += i.length();

        if(sequences_type == 0)
            this->kmer_size = (int)(log(this->genes_lenght) / log(this->alphabet.size()));
        else
            this->kmer_size = (int)(log(this->genes_lenght) / log(4));

        std::cout << "gene length: " << this->genes_lenght << " kmer size " << this->kmer_size << std::endl;

        this->check_kvalue_file(filename, sequences_type, output, logfile);
    }

    std::unordered_set<char> get_alphabet() {
        return this->alphabet;
    }

    [[nodiscard]] int get_kmer_size() const {
        return this->kmer_size;
    }

private:
    std::unordered_set<char> alphabet;
    int kmer_size;
    int genes_lenght;

    void check_kvalue_file(const char* filename, int sequences_type, std::string output, std::string logfile) const {
        if(kvalue != this->kmer_size*2) {

            std::ofstream ofs("../include/kvalue.h", std::ofstream::trunc);
            ofs << "#define kvalue " << this->kmer_size*2;
            ofs.close();

            std::stringstream ss;
            ss << "g++ -O3 ../main.cpp && a.out" << " -f " << filename << " -s " << sequences_type << " -o " << output << " -l " << logfile;

            std::string command = ss.str();

            std::cout << command.c_str() << std::endl;
            //system(command.c_str());
            system("echo \"ciao\"");
            exit(0);
        }
    }


};

#endif //PANDELOS_PLUSPLUS_KVALUE_H
