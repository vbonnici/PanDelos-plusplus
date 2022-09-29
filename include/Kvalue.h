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

/*
 * Class that is responsible for calculating the correct k-mer size value on the basis of the alphabet found by the sequences
 */
class Kvalue {
public:

    /*
     * The constructor calculates the alphabet, determines the value of k-mer size and checks that the configuration file is adequate for it
     *
     * @param[in] const std::vector<std::string>&
     * @param[in] const char*
     * @param[in] const int
     * @param[in] const std::string&
     * @param[in] const std::string&
     * @param[in] std::ofstream*
     */
    explicit Kvalue(const std::vector<std::string>& input_sequences, const char* filename, const int sequences_type, const std::string& output, const std::string& logfile, std::ofstream* log_stream) {
        if(log_stream->bad())
            throw std::runtime_error("the log stream's badbit error state flag is set");

        if(input_sequences.empty())
            throw std::runtime_error("input sequences vector is empty");

        this->log_stream = log_stream;
        this->input_sequences = &input_sequences;

        /*
         * For each gene sequence, it inserts each character of which it is composed into an unordered set
         * (if not yet present)
         */
        for (auto &i: *this->input_sequences)
            for(char a : i)
                this->alphabet.insert(a);

        *this->log_stream << "Identified alphabet: " << std::endl;

        for(auto &i : this->alphabet)
            *this->log_stream << i << " ";
        *this->log_stream << std::endl;

        this->genes_lenght = 0;
        for(auto &i : *this->input_sequences)
            this->genes_lenght += i.length();

        if(sequences_type == 0)
            this->kmer_size = (int)(log(this->genes_lenght) / log(this->alphabet.size()));
        else
            this->kmer_size = (int)(log(this->genes_lenght) / log(4));

        *this->log_stream << "Gene length: " << this->genes_lenght << std::endl;
        *this->log_stream << "Kmer size: " << this->kmer_size << std::endl;

        this->check_kvalue_file(filename, sequences_type, output, logfile);
    }

    /*
     * Getter which returns the calculated alphabet
     *
     * @param[out] std::unordered_set<char>
     */
    [[nodiscard]] std::unordered_set<char> get_alphabet() const {
        return this->alphabet;
    }

    /*
    * Getter which returns the k-mer size
    *
    * @param[out] int
    */
    [[nodiscard]] int get_kmer_size() const {
        return this->kmer_size;
    }

private:
    std::unordered_set<char> alphabet;
    int kmer_size;
    int genes_lenght;
    std::ofstream* log_stream;
    const std::vector<std::string>* input_sequences;

    /*
     * A software class uses std::bitset. This data structure requires a size of the bits to be stored which must be
     * known at compile-time.
     * To respect this constraint and guarantee a dynamic choice of the size of the bits to be stored, after calculating
     * the k-mer size, this method verifies that a preprocessing directive (of type define) is adequate for this value.
     *
     * Alternatively, the configuration file that contains this directive is overwritten with an appropriate value and the
     * program is recompiled and launched with the same original options chosen by the user
     *
     * @param[in] const char*
     * @param[in] const int
     * @param[in] const std::string&
     * @param[in] const std::string&
     */
    void check_kvalue_file(const char* filename, const int sequences_type, const std::string& output, const std::string& logfile) const {
        int kvalue_size_expected;

        if(sequences_type == 0)
            kvalue_size_expected = this->kmer_size*3*2; //2 bits for each nucleotide
        else
            kvalue_size_expected = this->kmer_size*2;


        if(kvalue != kvalue_size_expected) {
            std::ofstream ofs;

            std::filesystem::path cwd = std::filesystem::current_path() / "conf/conf.h";
            ofs = std::ofstream(cwd.string(), std::ofstream::trunc);

            if (!ofs) {
                *this->log_stream << "file opening error" << std::endl;
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