#ifndef FUNCIONES_HPP
#define FUNCIONES_HPP

#include <string>
#include <utility> 
#include <string>

std::string find_longest_palindrome(const std::string& s);
std::pair<std::string, std::string> read_fasta(const std::string& filename);

std::string codon_transformer(const std::string& seq);

#endif 
