#ifndef PANDELOS_PLUSPLUS_PANGENEIDATA_H
#define PANDELOS_PLUSPLUS_PANGENEIDATA_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <cstring>
#include <iterator>
#include "lib/Helper.h"

/*
 * .faa file reader in standard format
 *
 * The class is necessarily instantiated with a string that represents the directory of the .faa file to read.
 * The file is stored in a buffer and the data structures are built containing the components of the file divided by type
 * @see README.me
 * The class offers getters to get the single data structures and helper methods to print them on standard output
 *
 */
class PangeneIData {
public:
    //TODO: completa documentazione del costruttore dopo aver ottimizzato le strutture dati necessarie
    /*
     * Opens the file and stores in data structures
     *
     * @param[in] const char* filename
     */
    explicit PangeneIData(const char* filename, std::ofstream* log_stream) {
        this->log_stream = log_stream;
        this->open_file(filename);

        bool nameline = true;
        std::string genome_name, sequence_name, sequence_description;
        std::string genome_name_old;
        std::vector<int> temp_vector;
        int counter_sequence = 0;
        bool new_genome = false;

        while (this->is_valid()) {
            if (nameline) {
                genome_name = this->next_string();
                sequence_name = this->next_string();
                sequence_description = this->next_string();
            } else {
                std::string sequence = this->next_string();
                this->sequences.emplace_back(sequence);
                this->sequences_name.emplace_back(sequence_name);
                this->sequences_description.push_back(sequence_description);

                if(genome_name != genome_name_old) {
                    if(new_genome) {
                        this->genome_sequencesid.emplace_back(temp_vector);
                        temp_vector.clear();
                    }

                    temp_vector.emplace_back(counter_sequence);
                    genome_name_old = genome_name;
                    new_genome = true;
                } else {
                    temp_vector.emplace_back(counter_sequence);
                }
                ++counter_sequence;
            }

            nameline = !nameline;
        }

        this->genome_sequencesid.emplace_back(temp_vector);

        this->genes_id_interval.reserve(this->genome_sequencesid.size());

        for(int i = 0; i < this->genome_sequencesid.size(); ++i) {
            auto genes = this->genome_sequencesid.operator[](i);

            auto min = min_element(std::begin(genes), std::end(genes));
            auto max = max_element(std::begin(genes), std::end(genes));

            this->genes_id_interval.emplace_back(std::make_pair(*min, *max));
        }

        *this->log_stream << "File di input letto correttamente" << std::endl;
    }

    /*
     * Kvalue that takes care of printing the vector containing the names of the genomes
     */
    void print_genomes_names() {
        Helper::simple_container_print<std::string, std::vector<std::string>::iterator>
            (*this->log_stream, this->genomes_names.begin(), this->genomes_names.end(), "\n");
    }

    /*
     * Kvalue that takes care of printing the vector containing the sequences description
     */
    void print_sequences_description() {
        Helper::simple_container_print<std::string, std::vector<std::string>::iterator>
                (*this->log_stream, this->sequences_description.begin(), this->sequences_description.end(), "\n");
    }

    /*
     * Kvalue that takes care of printing the vector containing the sequences
     */
    void print_sequences() {
        Helper::simple_container_print<std::string, std::vector<std::string>::iterator>
                (*this->log_stream, this->sequences.begin(), this->sequences.end(), "\n");
    }

    /*
     * Getter which returns the vector containing the names of all genomes contained in the .faa file
     *
     * @param[out] std::vector<std::string>&
     */
    std::vector<std::string>& get_genomes_names() {
        return this->genomes_names;
    }

    /*
     * Getter which returns the vector containing the names of all genes contained in the .faa file
     *
     * @param[out] std::vector<std::string>&
     */
    std::vector<std::string>& get_sequences_name() {
        return this->sequences_name;
    }

    /*
     * Getter which returns the vector containing the decription of all genes contained in the .faa file
     *
     * @param[out] std::vector<std::string>&
     */
    std::vector<std::string>& get_sequences_description() {
        return this->sequences_description;
    }

    /*
     * Getter which returns the vector containing the sequences of all genes contained in the .faa file
     *
     * @param[out] std::vector<std::string>&
     */
    std::vector<std::string>& get_sequences() {
        return this->sequences;
    }

    std::vector<std::vector<int>>& get_genome_sequencesid() {
        return this->genome_sequencesid;
    }

    std::vector<std::pair<int, int>>& get_genes_id_interval() {
        return this->genes_id_interval;
    }

    /*
     * Method that encapsulates the call to fclose () function which dissociates the named stream from its underlying
     * and free () which frees the memory occupied by the buffer
     *
     * @note This method must be called explicitly from instances of the class
     *
     * @throws std::runtime_error if the stream resource cannot be released
     *
     * @throws std::runtime_error if there are errors in the memory deallocation
     */
    void close() {
        try {
            if(fclose(this->input_file) != 0)
                throw std::runtime_error(std::strerror(errno));

            free(this->buffer);

            if(this->buffer == NULL)
                throw std::runtime_error(std::strerror(errno));

        } catch(std::exception const& e) {
            *this->log_stream << "Exception: " << e.what() << "\n";
            exit(11);
        }
    }

private:
    FILE* input_file{};
    long file_size{};
    char* buffer{};
    long cursor{};

    std::ofstream* log_stream;
    std::vector<std::string> sequences;
    std::vector<std::string> sequences_name;
    std::vector<std::vector<int>> genome_sequencesid;
    std::vector<std::string> sequences_description;
    std::vector<std::string> genomes_names;
    std::vector<std::pair<int, int>> genes_id_interval;

    /*
     * This method takes care of opening the file, allocating the buffer and storing all the contents of the file in it
     *
     * @param[in] const char* filename
     *
     * @throws std::runtime_error if the file cannot be opened
     *
     * @throws std::runtime_error if you are unable to allocate the buffer
     *
     * @throws std::runtime_error if the size of the data read does not match the size of the file
     */
    void open_file(const char* filename) {
        size_t result;
        this->cursor = 0;

        try {
            this->input_file = fopen(filename, "rb");
            if (this->input_file == NULL)
                throw std::runtime_error(std::strerror(errno));

            // obtain file size:
            fseek(this->input_file, 0, SEEK_END);
            this->file_size = ftell(this->input_file);
            rewind(this->input_file);

            // allocate memory to contain the whole file:
            this->buffer = (char *) malloc(sizeof(char) * this->file_size);
            if (this->buffer == NULL)
                throw std::runtime_error(std::strerror(errno));

            // copy the file into the buffer:
            result = fread(this->buffer, 1, this->file_size, this->input_file);
            if (result != this->file_size)
                throw std::runtime_error("size fread != filesize");

        } catch(std::exception const& e) {
            *this->log_stream << "Exception: " << e.what() << "\n";
            exit(11);
        }

        //*this->log_stream << buffer << "\n";
    }

    /*
     * @param[out] true if the cursor is less than the length of the file
     * @param[out] false otherwise
     */
    [[nodiscard]] bool is_valid() const {
        return this->cursor < this->file_size;
    }

    /*
     * Each time this method is called, a "string" is extracted from the buffer. A string, in this case, is a contiguous
     * sequence of characters up to \ t, \ r, or \ n
     *
     * @param[out] char*
     */
    char* next_string() {

        while ((this->is_valid()) && ((this->buffer[this->cursor] == '\n') || (this->buffer[this->cursor] == '\t') || (this->buffer[this->cursor] == '\r'))) {
            ++this->cursor;
        }

        long ci = this->cursor;

        while ((this->is_valid()) && (this->buffer[this->cursor] != '\n') && (this->buffer[this->cursor] != '\t') && (this->buffer[this->cursor] != '\r')) {
            ++this->cursor;
        }

        if (this->is_valid()) {
            this->buffer[this->cursor] = '\0';
        }
        ++this->cursor;

        return this->buffer + ci;
    }
};


#endif //PANDELOS_PLUSPLUS_PANGENEIDATA_H