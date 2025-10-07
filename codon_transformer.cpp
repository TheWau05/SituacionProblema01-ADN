//tenemos_que_hacer_que_al leer tttaaac cambie de reading frame para el siguiente codon
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
    };

    std::string protein;
    for (size_t i = 0; i + 2 < seq.length(); i += 3) {
        std::string codon = seq.substr(i, 3);
        if (codon_table.find(codon) != codon_table.end()) {
            char aa = codon_table[codon];
            if (aa == '*') break;  // Stop codon
            protein += aa;
        }
    }
    return protein;
}