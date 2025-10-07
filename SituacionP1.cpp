#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "Funciones.hpp"

using namespace std;

/*
 * ALGORITMO: Busqueda ingenua de patron (Naive Pattern Search)
 * Complejidad: O(n*m) donde n = texto, m = patron
 * Uso: Encontrar todas las apariciones de un gen en el genoma completo
 * 
 * Nota: Se usa string::find() de STL que implementa Boyer-Moore o similar
 * optimizado, pero conceptualmente es busqueda de subcadena.
 */
vector<int> find_all_occurrences(const string& text, const string& pattern) {
    vector<int> positions;
    size_t pos = text.find(pattern, 0);
    while(pos != string::npos) {
        positions.push_back(pos);
        pos = text.find(pattern, pos + 1);
    }
    return positions;
}

// Funcion para obtener el codon desde una posicion en el genoma
string get_codon_at_position(const string& genome, int pos) {
    if (pos + 2 < genome.length()) {
        return genome.substr(pos, 3);
    }
    return "";
}

int main() {
    cout << "=== ANALISIS DEL GENOMA SARS-CoV-2 ===\n\n";
    
    // ========== PARTE 1: Encontrar indices de genes en el genoma ==========
    cout << "PARTE 1: Indices de genes en el genoma\n";
    cout << "========================================\n\n";
    
    auto [header_genome, genome] = read_fasta("SARS-COV-2-MN908947.3.txt");
    
    vector<pair<string, string>> genes = {
        {"gen-S.txt", "Gen S (Spike)"},
        {"gen-M.txt", "Gen M (Membrane)"},
        {"gen-ORF1AB.txt", "Gen ORF1AB"}
    };
    
    /*
     * ALGORITMO: Busqueda de subcadena
     * Para cada gen, busca su secuencia completa dentro del genoma
     */
    for (const auto& [file, name] : genes) {
        auto [header, gene_seq] = read_fasta(file);
        
        if (gene_seq.empty()) {
            cout << name << ": No se pudo leer\n\n";
            continue;
        }
        
        vector<int> positions = find_all_occurrences(genome, gene_seq);
        
        cout << name << ":\n";
        cout << "  Archivo: " << file << "\n";
        
        if (positions.empty()) {
            cout << "  No encontrado en el genoma\n";
        } else {
            cout << "  Indices de aparicion: ";
            for (size_t i = 0; i < positions.size(); i++) {
                cout << positions[i];
                if (i < positions.size() - 1) cout << ", ";
            }
            cout << "\n";
        }
        
        string first12 = gene_seq.substr(0, min<size_t>(12, gene_seq.size()));
        cout << "  Primeros 12 caracteres: " << first12 << "\n\n";
    }
    
    // ========== PARTE 2: Palindromos mas largos ==========
    cout << "\nPARTE 2: Palindromos mas largos en cada gen\n";
    cout << "=============================================\n\n";
    
    /*
     * ALGORITMO DE MANACHER
     * Implementado en find_longest_palindrome() en Funciones.cpp
     * Complejidad: O(n)
     * Encuentra el palindromo mas largo en una cadena
     * 
     * Funcionamiento:
     * 1. Transforma la cadena agregando caracteres especiales para 
     *    manejar palindromos de longitud par e impar uniformemente
     * 2. Usa el concepto de "centro" y "radio" para expandir palindromos
     * 3. Aprovecha la simetria de palindromos previamente encontrados
     *    para evitar comparaciones redundantes
     */
    ofstream palindrome_file("palindromos.txt");
    palindrome_file << "PALINDROMOS MAS LARGOS EN GENES SARS-CoV-2\n";
    palindrome_file << "==========================================\n\n";
    
    for (const auto& [file, name] : genes) {
        auto [header, gene_seq] = read_fasta(file);
        
        if (gene_seq.empty()) continue;
        
        string palindrome = find_longest_palindrome(gene_seq);
        
        cout << name << ":\n";
        cout << "  Longitud del palindromo mas largo: " << palindrome.length() << " nt\n";
        cout << "  Secuencia: " << palindrome << "\n\n";
        
        palindrome_file << name << ":\n";
        palindrome_file << "  Longitud: " << palindrome.length() << " nt\n";
        palindrome_file << "  Secuencia: " << palindrome << "\n\n";
    }
    
    palindrome_file.close();
    cout << "  [Palindromos guardados en palindromos.txt]\n\n";
    
    // ========== PARTE 3: Localizacion de proteinas en el genoma ==========
    cout << "\nPARTE 3: Localizacion de proteinas en el genoma\n";
    cout << "================================================\n\n";
    
    /*
     * ALGORITMO: Busqueda exhaustiva con traduccion en tiempo real
     * Complejidad: O(3 * n * m) donde n = genoma, m = proteina
     * 
     * Para cada proteina:
     * 1. Itera sobre las 3 fases de lectura posibles (frame 0, 1, 2)
     * 2. En cada posicion, traduce codones a aminoacidos
     * 3. Compara la traduccion con la secuencia de la proteina
     * 4. Si hay coincidencia completa, reporta la ubicacion
     * 
     * Optimizacion posible: Usar KMP o Boyer-Moore en el espacio de proteinas
     */
    vector<ProteinInfo> proteins = read_protein_into_file("seq-proteins.txt");
    cout << "Total de proteinas leidas: " << proteins.size() << "\n\n";
    
    map<string, char> codon_table = {
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
        {"TGC", 'C'}, {"TGT", 'C'}, {"TGA", '*'}, {"TGG", 'W'}
    };
    
    for (const auto& protein : proteins) {
        if (protein.sequence.length() < 4) continue;
        
        cout << "Proteina: " << protein.name << "\n";
        
        // Buscar en las 3 fases de lectura (reading frames)
        bool found = false;
        for (int frame = 0; frame < 3 && !found; frame++) {
            for (size_t i = frame; i + protein.sequence.length() * 3 < genome.length(); i += 3) {
                string translated = "";
                bool match = true;
                
                // Traducir segmento del genoma codon por codon
                for (size_t j = 0; j < protein.sequence.length() && match; j++) {
                    if (i + j * 3 + 2 >= genome.length()) {
                        match = false;
                        break;
                    }
                    
                    string codon = genome.substr(i + j * 3, 3);
                    if (codon_table.find(codon) != codon_table.end()) {
                        char aa = codon_table[codon];
                        if (aa == '*') break;  // Codon de paro
                        if (aa != protein.sequence[j]) {
                            match = false;
                        }
                    } else {
                        match = false;
                    }
                }
                
                if (match) {
                    found = true;
                    cout << "  Indices en genoma: " << i << " - " << (i + protein.sequence.length() * 3 - 1) << "\n";
                    cout << "  Fase de lectura: " << frame << "\n";
                    
                    // Primeros 4 aminoacidos
                    string first4aa = protein.sequence.substr(0, min<size_t>(4, protein.sequence.length()));
                    cout << "  Primeros 4 aminoacidos: " << first4aa << "\n";
                    
                    // Codones asociados
                    cout << "  Codones asociados: ";
                    for (size_t k = 0; k < 4 && k < protein.sequence.length(); k++) {
                        if (i + k * 3 + 2 < genome.length()) {
                            cout << genome.substr(i + k * 3, 3);
                            if (k < 3 && k < protein.sequence.length() - 1) cout << " ";
                        }
                    }
                    cout << "\n\n";
                    break;
                }
            }
        }
        
        if (!found) {
            cout << "  No encontrada en el genoma\n\n";
        }
    }
    
    // ========== PARTE 4: Comparacion Wuhan vs Texas ==========
    cout << "\nPARTE 4: Comparacion de genomas (Wuhan 2019 vs Texas 2020)\n";
    cout << "===========================================================\n\n";
    
    /*
     * ALGORITMO: Comparacion directa nucleotido por nucleotido
     * Complejidad: O(n) donde n = longitud del genoma
     * 
     * Para cada posicion:
     * 1. Compara nucleotidos en ambos genomas
     * 2. Si difieren, identifica el codon afectado
     * 3. Traduce ambos codones a aminoacidos
     * 4. Determina si es mutacion sinonima o no sinonima
     * 
     * Tipos de mutaciones:
     * - Sinonima: Cambia nucleotido pero produce el mismo aminoacido
     * - No sinonima: Cambia nucleotido y produce diferente aminoacido
     */
    
    auto [header_wuhan, genome_wuhan] = read_fasta("SARS-COV-2-MN908947.3.txt");
    auto [header_texas, genome_texas] = read_fasta("SARS-COV-2-MT106054.1.txt");
    
    cout << "Longitud genoma Wuhan: " << genome_wuhan.length() << " nt\n";
    cout << "Longitud genoma Texas: " << genome_texas.length() << " nt\n\n";
    
    size_t min_len = min(genome_wuhan.length(), genome_texas.length());
    vector<int> differences;
    
    // Encontrar diferencias
    for (size_t i = 0; i < min_len; i++) {
        if (genome_wuhan[i] != genome_texas[i]) {
            differences.push_back(i);
        }
    }
    
    cout << "Total de diferencias encontradas: " << differences.size() << "\n\n";
    
    if (differences.size() > 0) {
        cout << "Diferencias detalladas:\n";
        cout << "-----------------------\n\n";
        
        for (int pos : differences) {
            cout << "Posicion " << pos << ":\n";
            cout << "  Wuhan: " << genome_wuhan[pos] << "\n";
            cout << "  Texas: " << genome_texas[pos] << "\n";
            
            // Determinar posicion en codon (cual de los 3 nucleotidos es)
            int codon_start = (pos / 3) * 3;
            if (codon_start + 2 < min_len) {
                string codon_wuhan = genome_wuhan.substr(codon_start, 3);
                string codon_texas = genome_texas.substr(codon_start, 3);
                
                cout << "  Codon Wuhan: " << codon_wuhan;
                if (codon_table.find(codon_wuhan) != codon_table.end()) {
                    cout << " -> " << codon_table[codon_wuhan];
                }
                cout << "\n";
                
                cout << "  Codon Texas: " << codon_texas;
                if (codon_table.find(codon_texas) != codon_table.end()) {
                    cout << " -> " << codon_table[codon_texas];
                }
                cout << "\n";
                
                // Determinar si cambia el aminoacido
                if (codon_table.find(codon_wuhan) != codon_table.end() && 
                    codon_table.find(codon_texas) != codon_table.end()) {
                    if (codon_table[codon_wuhan] != codon_table[codon_texas]) {
                        cout << "  *** MUTACION NO SINONIMA: Cambia aminoacido ***\n";
                    } else {
                        cout << "  (Mutacion sinonima: mismo aminoacido)\n";
                    }
                }
            }
            cout << "\n";
        }
    }
    
    cout << "\n=== ANALISIS COMPLETADO ===\n";
    
    return 0;
}