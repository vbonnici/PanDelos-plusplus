#ifndef PANDELOS_PLUSPLUS_KVALUE_H
#define PANDELOS_PLUSPLUS_KVALUE_H

#include <ostream>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include "../include/kvalue.h"
#include <stdlib.h>
#include <filesystem>
#include "../include/global_options.h"

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

        this->sanitize_alphabet();

        *this->log_stream << "alfabeto calcolato: " << std::endl;

        //print alphabet
        for(auto &i : this->alphabet)
            *this->log_stream << i << " ";
        *this->log_stream << std::endl;

        this->genes_lenght = 0;
        for(auto &i : *sequences_input)
            this->genes_lenght += i.length();

        if(sequences_type == 0) {
            this->kmer_size = (int)(log(this->genes_lenght) / log(this->alphabet.size()));

            if(this->kmer_size % 3 != 0)
                this->kmer_size = this->roundUp(this->kmer_size, 3);
        }
        else
            this->kmer_size = (int)(log(this->genes_lenght) / log(4));

        *this->log_stream << "1 - gene length: " << this->genes_lenght << " kmer size " << this->kmer_size << std::endl;

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

        if(sequences_type == 1)
            kvalue_size_expected = this->kmer_size*2;   //2 bit per ogni nucleotide
        else
            kvalue_size_expected = (this->kmer_size/3)*2; //triplette e 2 bit per ogni nucleotide

        if(kvalue != kvalue_size_expected) {
            std::ofstream ofs;

            if(debug)
                ofs = std::ofstream("../include/kvalue.h", std::ofstream::trunc);
            else {
                std::filesystem::path cwd = std::filesystem::current_path() / "include/kvalue.h";
                ofs = std::ofstream(cwd.string(), std::ofstream::trunc);
            }

            if (!ofs) {
                *this->log_stream << "errore apertura file" << std::endl;
                exit(11);
            }

            ofs << "#define kvalue " << kvalue_size_expected;
            ofs.close();


            std::stringstream ss;

            if(debug)
                ss << "g++ -w -std=c++17 -fopenmp -O3 ../main.cpp && ./a.out" << " -f " << filename << " -s " << sequences_type << " -o " << output << " -l " << logfile;
            else
                ss << "g++ -w -std=c++17 -fopenmp -O3 main.cpp && ./a.out" << " -f " << filename << " -s " << sequences_type << " -o " << output << " -l " << logfile;

            std::string command = ss.str();

            *this->log_stream << command.c_str() << std::endl;
            system(command.c_str());
            exit(0);
        }
    }

     int roundUp(int numToRound, int multiple) {
        if (multiple == 0)
            return numToRound;

        int remainder = numToRound % multiple;
        if (remainder == 0)
            return numToRound;

        return numToRound + multiple - remainder;
    }

    void sanitize_alphabet() {
        for(auto &i : this->alphabet) {
            std::string temp = reinterpret_cast<const char *>(i);
            if(temp.find_first_not_of("FLIMVSPTAY*HQNKDECWRG") != std::string::npos)
                this->alphabet.erase(i);
        }
    }

};

#endif //PANDELOS_PLUSPLUS_KVALUE_H
