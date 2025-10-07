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

        vector<ProteinInfo> proteins = read_protein_into_file("seq-proteins.txt");
        cout << proteins.size() << " proteinas leidas de seq-proteins.txt\n";

        for (const auto& protein : proteins) {
            string comparison = protein_compare(codon, protein.sequence);
            cout << "    * Comparación con " << protein.name << ": " << comparison << "\n";
            //auto [pos, segment] = find_protein_in_genome(seq, protein.sequence);
            //if (pos != -1) {
            //    cout << "      - Proteina encontrada en genoma en posición " << pos << ": " << segment << "\n";
            //} else {
            //    cout << "      - Proteina no encontrada en el genoma.\n";
            //}
        }
    }

    return 0;
}
