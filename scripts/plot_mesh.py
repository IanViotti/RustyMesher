import pandas as pd
import matplotlib.pyplot as plt

def plot_generated_mesh_scatter(csv_filepath):
    print(f"Lendo o arquivo: {csv_filepath}...")
    
    # 1. Carregar os dados do CSV
    try:
        df = pd.read_csv(csv_filepath)
    except FileNotFoundError:
        print(f"Erro: Arquivo '{csv_filepath}' não encontrado.")
        return

    # 2. Criar a figura
    plt.figure(figsize=(10, 8))

    # 3. Separar os dados para destaque
    # j == 0 é a fronteira interna (o perfil biconvexo)
    aerofolio = df[df['j'] == 0]
    malha_externa = df[df['j'] != 0]

    # 4. Plotar os pontos da malha como scatter (pontos azuis pequenos)
    # 's' é o tamanho do ponto (size) e 'alpha' dá uma leve transparência
    plt.scatter(malha_externa['x'], malha_externa['y'], 
                color='blue', s=5, alpha=0.5, label='Nós da Malha Computacional')

    # 5. Plotar os pontos do aerofólio (pontos pretos maiores)
    plt.scatter(aerofolio['x'], aerofolio['y'], 
                color='black', s=15, alpha=0.5, label='Superfície do Aerofólio')

    # 6. Configurações visuais do gráfico
    plt.axis('equal') # Mantém a proporção física correta x/y
    plt.title('Distribuição de Nós - Malha Tipo "O"')
    plt.xlabel('Coordenada X')
    plt.ylabel('Coordenada Y')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)

    print("Gerando o gráfico scatter...")
    plt.show()

if __name__ == "__main__":
    # Especifique o caminho exato do seu CSV
    arquivo_malha = "job_files/biconvex_mesh/malha_gerada.csv" 
    plot_generated_mesh_scatter(arquivo_malha)