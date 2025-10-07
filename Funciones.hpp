#ifndef FUNCIONES_HPP
#define FUNCIONES_HPP

#include <string>
#include <utility> 
#include <string>
#include <vector>

struct ProteinInfo {
    std::string name;
    std::string sequence;
};

std::string find_longest_palindrome(const std::string& s);
std::pair<std::string, std::string> read_fasta(const std::string& filename);
std::string codon_transformer(const std::string& seq);
int apariciones(const std::string &text, const std::string &sub);

std::string protein_compare(const std::string& seq1, const std::string& seq2);
std::vector<ProteinInfo> read_protein_into_file(const std::string& filename);
std::pair<int, std::string> find_protein_in_genome(const std::string& genome, const std::string& protein, int window_size = 3000);

#endif 
