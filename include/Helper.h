#ifndef PANDELOS_PLUSPLUS_HELPER_H
#define PANDELOS_PLUSPLUS_HELPER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <algorithm>

/*
 * Namespace that contains software utility functions and template functions
 */
namespace Helper {

    /*
     * Templatic function that allows to combine two maps into a single map which have the key of the same type.
     *
     * Preserving the keys and values of both
     *
     * @tparam KeyType type of key
     * @tparam LeftValue type of value of first map
     * @tparam RightValue type of value of second map
     *
     * @param[in] const std::map<KeyType, LeftValue>&
     * @param[in] const std::map<KeyType, RightValue>&
     *
     * @param[out] std::map<KeyType, std::pair<LeftValue, RightValue>>
     */
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

    /*
     * Comparison function necessary for the std::map data structure to store new bitsets
     */
    template<size_t sz> struct bitset_comparer {
        bool operator() (const std::bitset<sz> &b1, const std::bitset<sz> &b2) const {
            return b1.to_ulong() < b2.to_ulong();
        }
    };

    /*
    * Template function that allows to print on std::ostream a container identified by a pair of iterators
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

    /*
    * Template function that allows to print on std::ostream an unordered map
    * (start and end).
    *
    * @tparam map_key type of key
    * @tparam map_val type of value of map
    *
    * @param[in] std::ostream&
    * @param[in] const std::unordered_map<map_key, map_val>&
    * @param[in] const std::string&
    */
    template<typename map_key, typename map_val>
    void simple_unordered_map_print(std::ostream& ostr, const std::unordered_map<map_key, map_val>& _map, const std::string& delimiter) {
        for (const auto& item : _map)
            ostr << item.first << delimiter.c_str() << item.second << std::endl;
    }

    /*
    * Template function that allows to print on std::ostream an unordered map containing another unordered map
    * (start and end).
    *
    * @tparam map_key type of key
    * @tparam nested_map_key type of key of nested map
    * @tparam nested_map_value type of value of nested map
    *
    * @param[in] std::ostream&
    * @param[in] const std::unordered_map<map_key, std::unordered_map<nested_map_key, nested_map_value>>&
    * @param[in] const std::string&
    */
    template<typename map_key, typename nested_map_key, typename nested_map_value>
    void nested_map_print(std::ostream& ostr, const std::unordered_map<map_key, std::unordered_map<nested_map_key, nested_map_value>>& _map, const std::string& delimiter) {
        for (const auto& item : _map) {
            ostr << item.first << delimiter.c_str();

            simple_unordered_map_print<int, double>(ostr, item.second, delimiter);
        }
    }

    /*
    * Template function that allows to print on std::ostream an unordered map containing a pair
    * (start and end).
    *
    * @tparam map_key type of key
    * @tparam pair_val1 type of first item of pair
    * @tparam pair_val2 type of second item of pair
    *
    * @param[in] std::ostream&
    * @param[in] const std::unordered_map<map_key, std::pair<pair_val1, pair_val2>>&
    * @param[in] const std::string&
    */
    template<typename map_key, typename pair_val1, typename pair_val2>
    void pair_unordered_map_print(std::ostream& ostr, const std::unordered_map<map_key, std::pair<pair_val1, pair_val2>>& _map, const std::string& delimiter) {
        for (const auto& item : _map) {
            auto pair = item.second;
            ostr << item.first << delimiter.c_str() << pair.first << delimiter.c_str() << pair.second << std::endl;

        }
    }

    /*
    * Template function that allows to print on std::ostream, separately, the values of a pair
    * (start and end).
    *
    * @tparam type1 type of first item of pair
    * @tparam type2 type of second item of pair
    *
    * @param[in] std::ostream&
    * @param[in] const std::pair<type1, type2>
    * @param[in] const std::string&
    */
    template<typename type1, typename type2>
    void simple_pair_print(std::ostream& ostr, const std::pair<type1, type2> _pair, const std::string& delimiter) {
        ostr << std::get<0>(_pair) << delimiter.c_str() << std::get<1>(_pair) << std::endl;
    }

    /*
    * Template function that allows to print on std::ostream, separately, the first three values of a tuple
    * (start and end).
    *
    * @tparam type1 type of first item of tuple
    * @tparam type2 type of second item of tuple
    * @tparam type3 type of third item of tuple
    *
    * @param[in] std::ostream&
    * @param[in] const std::tuple<type1, type2, type3>
    * @param[in] const std::string&
    */
    template<typename type1, typename type2, typename type3>
    void simple_triple_print(std::ostream& ostr, const std::tuple<type1, type2, type3> _tuple, const std::string& delimiter) {
        ostr << std::get<0>(_tuple) << delimiter.c_str() << std::get<1>(_tuple) << delimiter.c_str() << std::get<2>(_tuple) << std::endl;
    }

    /*
     * Utility method that deals with transforming an amino acid into the corresponding nucleotide chain.
     * The transformation respects the official correspondence table between nucleotides and amino acids of prokaryotes
     */
    std::string aminoacid_to_nucleotides(const std::string& aminoacid) {
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