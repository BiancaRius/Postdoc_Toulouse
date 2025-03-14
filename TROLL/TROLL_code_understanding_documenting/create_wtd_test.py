def create_file_wtd(file_name="wtd_test.txt"):
    
    try:
        with open(file_name, 'w') as file:
            file.write("WTD_test\n")  # Escreve o cabeçalho
            for _ in range(5):
                file.write("1.5\n")  # Escreve o valor 1.5 em cada linha

        print(f"File '{file_name}' criado com sucesso.")
    except Exception as e:
        print(f"errot: {e}")

# Exemplo de uso
create_file_wtd()