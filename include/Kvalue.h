#ifndef PANDELOS_PLUSPLUS_KVALUE_H
#define PANDELOS_PLUSPLUS_KVALUE_H

#include <ostream>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include "../conf/conf.h"
#include <stdlib.h>
#include <filesystem>

class Kvalue {
public:

    /*
     * For each gene sequence, it inserts each character of which it is composed into an unordered set
     * (if not yet present)
     */
    explicit Kvalue(const std::vector<std::string>* sequences_input, const char* filename, int sequences_type, const std::string& output, const std::string& logfile, std::ofstream* log_stream) {
        this->log_stream = log_stream;

        for (auto &i: *sequences_input)
            for(char a : i)
                this->alphabet.insert(a);

        *this->log_stream << "Calculated alphabet: " << std::endl;

        for(auto &i : this->alphabet)
            *this->log_stream << i << " ";
        *this->log_stream << std::endl;

        this->genes_lenght = 0;
        for(auto &i : *sequences_input)
            this->genes_lenght += i.length();

        if(sequences_type == 0)
            this->kmer_size = (int)(log(this->genes_lenght) / log(this->alphabet.size()));
        else
            this->kmer_size = (int)(log(this->genes_lenght) / log(4));

        *this->log_stream << "Gene length: " << this->genes_lenght << std::endl;
        *this->log_stream << "Kmer size: " << this->kmer_size << std::endl;

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
    std::ofstream* log_stream;

    void check_kvalue_file(const char* filename, int sequences_type, const std::string& output, const std::string& logfile) const {
        int kvalue_size_expected;

        if(sequences_type == 0)
            kvalue_size_expected = this->kmer_size*3*2; //2 bit per ogni nucleotide
        else
            kvalue_size_expected = this->kmer_size*2;


        if(kvalue != kvalue_size_expected) {
            std::ofstream ofs;

            std::filesystem::path cwd = std::filesystem::current_path() / "conf/conf.h";
            ofs = std::ofstream(cwd.string(), std::ofstream::trunc);

            if (!ofs) {
                *this->log_stream << "errore apertura file" << std::endl;
                exit(11);
            }

            ofs << "#define kvalue " << kvalue_size_expected << std::endl;
            ofs.close();


            std::stringstream ss;

            ss << "g++ -w -std=c++17 -fopenmp -O3 src/cpp/main.cpp -o bin/pandelos_plus_plus.out && bin/pandelos_plus_plus.out" << " -f " << filename << " -s " << sequences_type << " -o " << output << " -l " << logfile;

            std::string command = ss.str();

            *this->log_stream << command.c_str() << std::endl;
            system(command.c_str());
            exit(0);
        }
    }
};

#endif //PANDELOS_PLUSPLUS_KVALUE_H
