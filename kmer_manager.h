#ifndef PANDELOS_PLUSPLUS_KMER_MANAGER_H
#define PANDELOS_PLUSPLUS_KMER_MANAGER_H

template<size_t sz> struct bitset_comparer {
    bool operator() (const std::bitset<sz> &b1, const std::bitset<sz> &b2) const {
        return b1.to_ulong() < b2.to_ulong();
    }
};


class kmer_manager {
public:
    void get_current_bitset() {
        if(!bitset_18.empty())
            this->get_bitset_18();
        else if (!bitset_20.empty())
            this->get_bitset_20();
    }

    void init_container(int kmer_size) {
        if(kmer_size == 9)
            this->init_bitset18();

        else if(kmer_size == 10)
            this->init_bitset20();
    }

    std::unordered_map<std::string, std::map<std::bitset<18>, int, bitset_comparer<18>>>& get_bitset_18() {
        return this->bitset_18;
    }

    std::unordered_map<std::string, std::map<std::bitset<20>, int, bitset_comparer<20>>>& get_bitset_20() {
        return this->bitset_20;
    }

    void add_kmer(std::string& sequence) {
        if(!bitset_18.empty())
            init_bitset18(sequence);
        else if (!bitset_20.empty())
            init_bitset20(sequence);
    }


private:
    std::unordered_map<std::string, std::map<std::bitset<18>, int, bitset_comparer<18>>> bitset_18;
    std::unordered_map<std::string, std::map<std::bitset<20>, int, bitset_comparer<20>>> bitset_20;


    void init_bitset18() {
        std::map<std::bitset<18>, int, bitset_comparer<18>> temp;
        std::bitset<18> bitset_temp;
        bitset_temp.set(false);
        temp.insert(std::make_pair(bitset_temp, 0));

        this->bitset_18.insert(std::make_pair("", temp));
    }

    void init_bitset18(std::string& sequence) {
        std::map<std::bitset<18>, int, bitset_comparer<18>> temp;
        std::bitset<18> bitset_temp;
        bitset_temp.set(false);
        temp.insert(std::make_pair(bitset_temp, 0));

        this->bitset_18.insert(std::make_pair(sequence, temp));
    }

    void init_bitset20() {
        std::map<std::bitset<20>, int, bitset_comparer<20>> temp;
        std::bitset<20> bitset_temp;
        bitset_temp.set(false);
        temp.insert(std::make_pair(bitset_temp, 0));

        this->bitset_20.insert(std::make_pair("", temp));
    }

    void init_bitset20(std::string& sequence) {
        std::map<std::bitset<20>, int, bitset_comparer<20>> temp;
        std::bitset<20> bitset_temp;
        bitset_temp.set(false);
        temp.insert(std::make_pair(bitset_temp, 0));

        this->bitset_20.insert(std::make_pair(sequence, temp));
    }

};

#endif //PANDELOS_PLUSPLUS_KMER_MANAGER_H
