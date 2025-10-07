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

std::string find_longest_palindrome(const std::string& s) {
    if (s.empty()) return "";

    
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

