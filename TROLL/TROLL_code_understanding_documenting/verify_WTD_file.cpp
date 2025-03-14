#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::string file_path = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/Paracou_input_pedology_WTD.txt";
    std::ifstream file(file_path); // Abre o arquivo para leitura

    if (!file) {
        std::cerr << "Erro ao abrir o arquivo: " << file_path << std::endl;
        return 1; // Retorna erro
    }

    std::string line;
    int line_count = 0;

    while (std::getline(file, line) && line_count < 10) { // Lê e imprime as primeiras 10 linhas
        std::cout << line << std::endl;
        line_count++;
    }

    file.close(); // Fecha o arquivo
    return 0;
}
