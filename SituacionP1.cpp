#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

//tenemos_que_hacer_que_al leer tttaaac cambie de reading frame para el siguiente codon
//cada que se detecte tttaaac cambia

pair<string, string> read_fasta(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error abriendo archivo: " << filename << endl;
        return {"", ""};
    }
    string header, line, sequence;
    getline(file, header);
    while (getline(file, line)) {
        for (char c : line) {
            if (!isspace((unsigned char)c)) {
                char up = toupper(c);
                if (up == 'U') up = 'T';
                sequence.push_back(up);
            }
        }
    }
    return {header, sequence};
}

string manacher_longest_palindrome(const string &s) {
    if (s.empty()) return "";
    string t = "^#";
    for (char c : s) { t.push_back(c); t.push_back('#'); }
    t.push_back('$');
    int n = t.size();
    vector<int> p(n, 0);
    int center = 0, right = 0;
    int maxLen = 0, centerIndex = 0;
    for (int i = 1; i < n - 1; i++) {
        int mirror = 2 * center - i;
        if (i < right) p[i] = min(right - i, p[mirror]);
        while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) ++p[i];
        if (i + p[i] > right) { center = i; right = i + p[i]; }
        if (p[i] > maxLen) { maxLen = p[i]; centerIndex = i; }
    }
    int start = (centerIndex - maxLen) / 2;
    return s.substr(start, maxLen);
}

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

        string pal = manacher_longest_palindrome(seq);
        cout << "  - Palíndromo más largo: " << (pal.empty() ? "(ninguno)" : pal) << "\n";
        cout << "  - Caracteres: " << pal.size() << " nt\n" << "\n";

       
    }
    return 0;
}
