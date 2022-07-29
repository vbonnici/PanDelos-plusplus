#ifndef PANDELOS_PLUSPLUS_HELPER_H
#define PANDELOS_PLUSPLUS_HELPER_H

namespace Helper {

    template <typename KeyType, typename LeftValue, typename RightValue>
    std::map<KeyType, std::pair<LeftValue, RightValue>>
    UnionMaps(const std::map<KeyType, LeftValue>& map1, const std::map<KeyType, RightValue>& map2) {

        std::map<KeyType, std::pair<LeftValue, RightValue>> result;
        auto it_map1  = map1.begin();
        auto it_map2 = map2.begin();

        while (it_map1 != map1.end() && it_map2 != map2.end()) {
            if (it_map1->first < it_map2->first) {
                result.insert(std::make_pair(it_map1->first, std::make_pair(it_map1->second, 0)));
                it_map1++;
            } else if (it_map1->first == it_map2->first) {
                result.insert(std::make_pair(it_map1->first, std::make_pair(it_map1->second, it_map2->second)));
                it_map1++;
                it_map2++;
            } else if (it_map1->first > it_map2->first) {
                result.insert(std::make_pair(it_map2->first, std::make_pair(0, it_map2->second)));
                it_map2++;
            }
        }

        while (it_map1 != map1.end()) {
            result.insert(std::make_pair(it_map1->first, std::make_pair(it_map1->second, 0)));
            it_map1++;
        }

        while (it_map2 != map2.end()) {
            result.insert(std::make_pair(it_map2->first, std::make_pair(0, it_map2->second)));
            it_map2++;
        }

        return result;
    }

    template<size_t sz> struct bitset_comparer {
        bool operator() (const std::bitset<sz> &b1, const std::bitset<sz> &b2) const {
            return b1.to_ulong() < b2.to_ulong();
        }
    };

    /*
    * Template function that allows you to print on std::ostream a container identified by a pair of iterators
    * (start and end).
    *
    * @tparam T type of data that the container stores
    * @tparam InputIterator
    *
    * @param[in] std::ostream& ostr
    * @param[in] InputIterator itbegin
    * @param[in] InputIterator itend
    * @param[in] const std::string& delimiter
    */
    template<typename T, typename InputIterator>
    void simple_container_print(std::ostream& ostr, InputIterator itbegin, InputIterator itend, const std::string& delimiter) {
        std::copy(itbegin, itend, std::ostream_iterator<T>(ostr, delimiter.c_str()));
    }

    template<typename map_key, typename map_val>
    void simple_unordered_map_print(std::ostream& ostr, const std::unordered_map<map_key, map_val>& _map, const std::string& delimiter) {
        for (const auto& item : _map)
            ostr << item.first << delimiter.c_str() << item.second << std::endl;
    }

    template<typename map_key, typename nested_map_key, typename nested_map_value>
    void nested_unordered_map_print(std::ostream& ostr, const std::unordered_map<map_key, std::unordered_map<nested_map_key, nested_map_value>>& _map, const std::string& delimiter) {
        for (const auto& item : _map) {
            ostr << item.first << delimiter.c_str();

            simple_unordered_map_print<std::string, int>(ostr, item.second, delimiter);
        }
    }

    template<typename map_key, typename pair_val1, typename pair_val2>
    void pair_unordered_map_print(std::ostream& ostr, const std::unordered_map<map_key, std::pair<pair_val1, pair_val2>>& _map, const std::string& delimiter) {
        for (const auto& item : _map) {
            auto pair = item.second;
            ostr << item.first << delimiter.c_str() << pair.first << delimiter.c_str() << pair.second << std::endl;

        }
    }

    template<typename type1, typename type2, typename type3>
    void simple_triple_print(std::ostream& ostr, const std::tuple<type1, type2, type3> _tuple, const std::string& delimiter) {
        ostr << std::get<0>(_tuple) << delimiter.c_str() << std::get<1>(_tuple) << delimiter.c_str() << std::get<2>(_tuple) << std::endl;
    }

    std::string aminoacid_to_nucleotides(std::string& aminoacid) {
        std::string kmer;

        for(int a = 0; a < aminoacid.length(); ++a) {

            if(reinterpret_cast<char>(aminoacid[a]) == 'F')
                kmer += "TTT";

            else if(reinterpret_cast<char>(aminoacid[a]) == 'L')
                kmer += "TTA";

            else if(reinterpret_cast<char>(aminoacid[a]) == 'I')
                kmer += "ATT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'M')
                kmer += "ATG";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'V')
                kmer += "GTT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'S')
                kmer += "TCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'P')
                kmer += "CCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'T')
                kmer += "ACT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'A')
                kmer += "GCT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'Y')
                kmer += "TAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == '*')
                kmer += "TAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'H')
                kmer += "CAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'Q')
                kmer += "CAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'N')
                kmer += "AAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'K')
                kmer += "AAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'D')
                kmer += "GAT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'E')
                kmer += "GAA";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'C')
                kmer += "TGT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'W')
                kmer += "TGG";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'R')
                kmer += "CGT";

            else if (reinterpret_cast<char>(aminoacid[a]) == 'G')
                kmer += "GGG";
        }

        return kmer;
    }
}
#endif //PANDELOS_PLUSPLUS_HELPER_H
