#include "Funciones.hpp"
#include <iostream>
#include <fstream>
#include <cctype>

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
                sequence += up;  // ✅ Aquí agregamos cada carácter a la secuencia
            }
        }
    }

    return {header, sequence};
}
