#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "Funciones.hpp"

using namespace std;

int main() {

    vector<string> files = {
        "gen-S.txt",
        "gen-M.txt",
        "gen-ORF1AB.txt"
    };

    for (size_t i = 0; i < files.size(); i++) {
        auto [header, seq] = read_fasta(files[i]); 
        cout << "[" << (i+1) << "] " << files[i] << "\n";
        if (seq.empty()) { cout << "  > Secuencia vacía.\n\n"; continue; }

        string first12 = seq.substr(0, min<size_t>(12, seq.size()));
        cout << "  - Primeros 12 nt: " << first12 << "\n";

        string pal = find_longest_palindrome(seq);
        cout << "  - Palindromo mas largo: " << (pal.empty() ? "(ninguno)" : pal) << "\n";
        cout << "  - Caracteres: " << pal.size() << " nt\n" << "\n";

        string codon = codon_transformer(seq);
        cout << "  - Proteina traducida: " << (codon.empty() ? "(ninguna)" : codon) << "\n";
        cout << "  - Aminoacidos: " << codon.size() << " aa\n" << "\n";
    }

    return 0;
}
