#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream infile("Paracou_input_pedology_WTD.txt"); // Altere para o nome do seu arquivo

    // Verifica se o arquivo foi aberto com sucesso
    if (!infile) {
        std::cerr << "Erro ao abrir o arquivo." << std::endl;
        return 1;
    }

    std::string line;

    // Lê o arquivo linha por linha
    while (std::getline(infile, line)) {
        std::cout << line << std::endl; // Imprime cada linha no console
    }

    infile.close(); // Fecha o arquivo
    return 0;
}
