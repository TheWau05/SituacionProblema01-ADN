#include <iostream>
#include <string>
#include <vector>
using namespace std;

// Cuenta cuántas veces aparece 'sub' dentro de 'text'
int apariciones(const string &text, const string &sub) {
    if (sub.empty()) return 0;

    int count = 0;
    size_t pos = 0;

    while ((pos = text.find(sub, pos)) != string::npos) {
        count++;
        pos += sub.length(); 
    }

    return count;
}
