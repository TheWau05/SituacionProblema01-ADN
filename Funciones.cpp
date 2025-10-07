#include "Funciones.hpp"
#include <iostream>
#include <fstream>
#include <cctype>
#include <string>
#include <vector>
#include <map>

std::pair<std::string, std::string> read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {"", ""};
    }

    std::string header, line, sequence;
    
    // Leer el encabezado (línea que empieza con '>')
    getline(file, header);

    // Leer el resto de líneas (la secuencia)
    while (getline(file, line)) {
        for (char c : line) {
            if (!isspace(static_cast<unsigned char>(c))) {
                char up = toupper(static_cast<unsigned char>(c));
                if (up == 'U') up = 'T'; 
                sequence += up;  
            }
        }
    }

    return {header, sequence};
}

// See header file for function documentation.
std::string find_longest_palindrome(const std::string& s) {
    if (s.empty()) return "";

    // Transform the string to handle even-length palindromes easily.
    // e.g., "aba" -> "^#a#b#a#$"
    std::string t = "^#";
    for (char c : s) {
        t.push_back(c);
        t.push_back('#');
    }
    t.push_back('$');

    int n = t.size();
    std::vector<int> p(n, 0);
    int center = 0, right = 0;
    int maxLen = 0, centerIndex = 0;

    for (int i = 1; i < n - 1; i++) {
        int mirror = 2 * center - i;
        if (i < right) {
            p[i] = std::min(right - i, p[mirror]);
        }
        // Attempt to expand palindrome centered at i
        while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) {
            p[i]++;
        }
        // If palindrome centered at i expands past `right`,
        // adjust center and right boundary
        if (i + p[i] > right) {
            center = i;
            right = i + p[i];
        }
        // Update the longest palindrome found so far
        if (p[i] > maxLen) {
            maxLen = p[i];
            centerIndex = i;
        }
    }

    int start = (centerIndex - maxLen) / 2;
    return s.substr(start, maxLen);
}




using namespace std;

// Función para transformar una secuencia de ADN en una secuencia de aminoácidos
// IN: secuencia de ADN (string)
// OUT: secuencia de aminoácidos (string)
std::string codon_transformer(const std::string& seq) {
    std::map<std::string, char> codon_table = {
        {"ATA", 'I'}, {"ATC", 'I'}, {"ATT", 'I'}, {"ATG", 'M'},
        {"ACA", 'T'}, {"ACC", 'T'}, {"ACG", 'T'}, {"ACT", 'T'},
        {"AAC", 'N'}, {"AAT", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGC", 'S'}, {"AGT", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"CTA", 'L'}, {"CTC", 'L'}, {"CTG", 'L'}, {"CTT", 'L'},
        {"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCT", 'P'},
        {"CAC", 'H'}, {"CAT", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGA", 'R'}, {"CGC", 'R'}, {"CGG", 'R'}, {"CGT", 'R'},
        {"GTA", 'V'}, {"GTC", 'V'}, {"GTG", 'V'}, {"GTT", 'V'},
        {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCT", 'A'},
        {"GAC", 'D'}, {"GAT", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGT", 'G'},
        {"TCA", 'S'}, {"TCC", 'S'}, {"TCG", 'S'}, {"TCT", 'S'},
        {"TTC", 'F'}, {"TTT", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TAC", 'Y'}, {"TAT", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGC", 'C'}, {"TGT", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"TTTAAAC", 'X'}  // Custom codon to change reading frame
    };

    std::string protein;
    for (size_t i = 0; i + 2 < seq.length(); i += 3) {
        std::string codon = seq.substr(i, 3);
        if (codon_table.find(codon) != codon_table.end()) { // usa find para encontrar llave del hash
            char aa = codon_table[codon];
            if (aa == '*') break;  // Stop codon
            protein += aa;
        }
    }
    return protein;
}


std::string protein_compare(const std::string& seq1, const std::string& seq2) {
    size_t len1 = seq1.length();
    size_t len2 = seq2.length();
    size_t minLen = std::min(len1, len2);

    size_t matches = 0;
    for (size_t i = 0; i < minLen; ++i) {
        if (seq1[i] == seq2[i]) {
            ++matches;
        }
    }

    double similarity = (static_cast<double>(matches) / minLen) * 100.0;

    return "Similarity: " + std::to_string(similarity) + "% (" + std::to_string(matches) + " out of " + std::to_string(minLen) + " amino acids match)";
}




std::vector<ProteinInfo> read_protein_into_file(const std::string& filename) {
    std::vector<ProteinInfo> proteins;
    std::ifstream file(filename);
    std::string line;
    ProteinInfo current;
    
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!current.sequence.empty()) {
                proteins.push_back(current);
            }
            current.name = line.substr(1);
            current.sequence = "";
        } else {
            current.sequence += line;
        }
    }
    if (!current.sequence.empty()) {
        proteins.push_back(current);
    }
    return proteins;
}

std::pair<int, std::string> find_protein_in_genome(const std::string& genome, const std::string& protein, int window_size ) {
    for (size_t i = 0; i + 2 < genome.length(); i++) {
        // Try all three reading frames at this position
        for (int frame = 0; frame < 3; frame++) {
            if (i + frame + window_size > genome.length()) continue;
            
            std::string segment = genome.substr(i + frame, window_size);
            std::string translated = codon_transformer(segment);
            
            size_t found = translated.find(protein.substr(0, 50)); // Check first 50 amino acids
            if (found != std::string::npos) {
                // Verify complete match
                if (translated.substr(found).find(protein) == 0) {
                    return {i + frame, genome.substr(i + frame, protein.length() * 3)};
                }
            }
        }
    }
    return {-1, ""};
}
