#ifndef PANDELOS_PLUSPLUS_HELPER_H
#define PANDELOS_PLUSPLUS_HELPER_H
#include <iostream>
#include <string>
#include "argagg/argagg.hpp"

using argagg::parser_results;
using argagg::parser;
parser argparser {
        {
                {
                        "help",
                        {
                                "-h",
                                "--help"
                        },
                        "Print help menu",
                        0
                }, {
                        "filename",
                        {
                                "-f",
                                "--filename"
                        },
                        ".faa filename",
                        1
                }, {
                        "sequences_type",
                        {
                                "-s",
                                "--stype"
                        },
                        "0 amino acids, 1 nucleotides",
                        1
                }, {
                        "output",
                        {
                                "-o",
                                "--output"
                        },
                        ".net filename",
                        1
                }, {
                        "log",
                        {
                                "-l",
                                "--log"
                        },
                        ".txt filename",
                        1
                },
        }
};

class Helper {
public:
    void parse_arguments(int argc, char* argv[]) {
        this->do_args_parsing(argc, argv);
    }

    const char* get_filename() {
        return this->filename;
    }

    [[nodiscard]] int get_sequences_type() const {
        return this->sequences_type;
    }

    std::string get_output() {
        return this->output;
    }

    std::string get_log() {
        return this->log;
    }

private:
    const char* filename;
    int sequences_type;
    std::string output;
    std::string log;

    void do_args_parsing(int argc, char* argv[]) {

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << argv[0] << " 2.0\n" << '\n' << '\n' << argparser << '\n'
                << "Encountered exception while parsing arguments: " << e.what()
                << '\n';
            exit(1);
        }

        if (args["help"]) {
            std::cerr << argparser;
            exit(1);
        }

        if (!args["filename"]) {
            std::cerr << "'--filename' ('-f') option required\n";
            exit(1);
        } else {
            this->filename = args["filename"].as<const char*>();

            if(!(args["filename"].as<std::string>().substr(args["filename"].as<std::string>().find_last_of(".") + 1) == "faa")) {
                std::cerr << "the input file must have a .faa extension\n";
                exit(1);
            }
        }

        if (!args["output"]) {
            std::cerr << "'--output' ('-o') option required\n";
            exit(1);
        } else {
            this->output = args["output"].as<std::string>();

            if(!(this->output.substr(this->output.find_last_of(".") + 1) == "net")) {
                std::cerr << "the output file must have a .net extension\n";
                exit(1);
            }
        }

        if (!args["log"]) {
            std::cerr << "'--log' ('-l') option required\n";
            exit(1);
        } else {
            this->log = args["log"].as<std::string>();

            if(!(this->log.substr(this->log.find_last_of(".") + 1) == "txt")) {
                std::cerr << args.program << "the log file must have a .txt extension\n";
                exit(1);
            }
        }

        if (!args["sequences_type"]) {
            std::cerr << args.program << "'--sequences_type' ('-stype') option required\n";
            exit(1);
        } else {
            if(args["sequences_type"].as<int>() == 0 || args["sequences_type"].as<int>() == 1)
                this->sequences_type = args["sequences_type"].as<int>();
            else {
                std::cerr << "0 amino acids, 1 nucleotides\n";
                exit(1);
            }
        }
    }
};
#endif //PANDELOS_PLUSPLUS_HELPER_H
