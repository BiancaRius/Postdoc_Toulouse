def create_file_wtd(file_name="wtd_test.txt"):
    valores = [1.0, 1.5, 2.0, 2.5, 3.0] 
    try:
        with open(file_name, 'w') as file:
            file.write("WTD_test\n")  # Escreve o cabeçalho
            for valor in valores:
                file.write(f"{valor}\n")  

        print(f"File '{file_name}' criado com sucesso.")
    except Exception as e:
        print(f"errot: {e}")

# Exemplo de uso
create_file_wtd()