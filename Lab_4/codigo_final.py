from Lab_1.lab1 import hamming_simulation 
from Lab_2.Lab2 import LDPC_decode 
from Lab_3.test_bibs import generate_p_from_ebn0_db, run_simulation
import numpy as np
import time
import matplotlib.pyplot as plt
import os
from scipy.stats import norm




def plot_ber_comparison(eb_n0_db_axis, ber_data_dict, filename="ber_all_comparison.png"):
    """
    Plota múltiplas curvas de BER vs Eb/N0.

    Args:
        eb_n0_db_axis (list): Valores de Eb/N0 em dB para o eixo X.
        ber_data_dict (dict): Dicionário onde chaves são labels (ex: "LDPC (3,7) LLR")
                              e valores são as listas de BER correspondentes.
                              "BPSK Não Codificado (Teórico)" é uma chave especial.
        filename (str): Nome do arquivo para salvar o gráfico.
    """
    plt.figure(figsize=(12, 8))
    markers = ['o', '^', 'D', 'x', 's', 'v', '<', '>'] # Mais marcadores
    linestyles = ['-', ':', '-.', '--', '-', ':'] # Mais estilos
    
    # Usar um ciclo de cores para melhor distinção se houver muitas curvas
    # cmap = plt.get_cmap('tab10') # tab10 tem 10 cores distintas
    # colors = [cmap(i) for i in np.linspace(0, 1, len(ber_data_dict) + 1)] # +1 para BPSK
    
    idx = 0
    for label, ber_values in ber_data_dict.items():
        # current_color = colors[idx]
        if label == "BPSK Não Codificado (Teórico)":
            eb_n0_linear_uncoded = [10**(x/10.0) for x in eb_n0_db_axis if x is not None]
            actual_ber_uncoded = [norm.sf(np.sqrt(2 * ebn0_lin)) if ebn0_lin > 0 else 0.5 for ebn0_lin in eb_n0_linear_uncoded]
            # Certifique-se que actual_ber_uncoded tem o mesmo tamanho que eb_n0_db_axis (ou os eixos válidos)
            # Se eb_n0_db_axis pode ter Nones, precisamos filtrar
            valid_ebn0_axis_for_bpsk = [x for x in eb_n0_db_axis if x is not None]
            if len(actual_ber_uncoded) == len(valid_ebn0_axis_for_bpsk):
                 plt.semilogy(valid_ebn0_axis_for_bpsk, actual_ber_uncoded, marker='s', linestyle='--', label=label) # color=current_color)
            else:
                print(f"Aviso: Não foi possível plotar BPSK teórico devido a incompatibilidade de tamanho de eixo/dados.")

        elif ber_values is not None and len(ber_values) == len(eb_n0_db_axis):
            plt.semilogy(eb_n0_db_axis, ber_values,
                         marker=markers[idx % len(markers)],
                         linestyle=linestyles[idx % len(linestyles)],
                         label=label) # color=current_color)
            idx += 1
        elif ber_values is not None:
            print(f"Aviso: Comprimento dos dados BER para '{label}' ({len(ber_values)}) não corresponde ao eixo Eb/N0 ({len(eb_n0_db_axis)}). Não será plotado.")
        # else: ber_values is None, não faz nada (como para BPSK placeholder inicial)


    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("Bit Error Rate (BER)")
    plt.title("Comparação de Desempenho de Códigos")
    plt.grid(True, which="both", ls="--")
    plt.legend(loc='best') # Tenta encontrar a melhor localização para a legenda
    plt.ylim(1e-6, 1.1e0) # Ajustado ylim max para ver melhor o 0.5
    if eb_n0_db_axis and any(val is not None for val in eb_n0_db_axis):
        valid_ebn0 = [x for x in eb_n0_db_axis if x is not None]
        if valid_ebn0:
             plt.xlim(min(valid_ebn0), max(valid_ebn0))
    else:
        plt.xlim(0, 5)

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError: # Se __file__ não estiver definido (ex: REPL)
        script_dir = os.getcwd()
    full_path = os.path.join(script_dir, filename)
    plt.savefig(full_path)
    print(f"Gráfico comparativo salvo em: {full_path}")
    # plt.show() # Comente para não bloquear em scripts, descomente para ver interativamente
    plt.close() # Fecha a figura da memória




# --- Configurações da Simulação ---
def main():
    # from Lab_2.Lab2 import LDPC_decode 
    # from Lab_1.lab1 import hamming_simulation 
    print('oi')
    # Parâmetros Comuns
    EB_N0_DB_RANGE = np.arange(0.0, 5.1, 0.25).tolist() # Passo de 0.5 para testes
    AMPLITUDE_BPSK = 1.0

    # --- Dicionário para armazenar todos os resultados BER ---
    all_ber_results_for_plot = {
        "BPSK Não Codificado (Teórico)": None # Será calculado e plotado pela função de plot
    }

    # --- 1. Gerar valores de p para simulações BSC ---
    print("Gerando probabilidades de erro de bit (p) para simulações BSC a partir de Eb/N0_dB...")
    channel_p_values = generate_p_from_ebn0_db(EB_N0_DB_RANGE)
    print(f"Valores de p gerados: {[f'{p_val:.3e}' for p_val in channel_p_values]}")
    print("-" * 50)

    # --- 2. Simulação Hamming (7,4) usando BSC com 'p' derivados ---
    print("INICIANDO SIMULAÇÃO HAMMING (7,4)")
    start_time_hamming = time.time()
    try:
        ber_hamming = hamming_simulation(vector_p=channel_p_values)
        all_ber_results_for_plot["Hamming (7,4) (BSC)"] = ber_hamming
        print(f"BER Hamming: {[f'{b:.3e}' for b in ber_hamming]}")
    except Exception as e:
        print(f"Erro durante a simulação Hamming: {e}")
    end_time_hamming = time.time()
    print(f"Simulação Hamming concluída em {end_time_hamming - start_time_hamming:.2f} segundos.")
    print("-" * 50)

    # --- 3. Simulação LDPC com Bit-Flipping usando BSC com 'p' derivados ---
    DV_LDPC_BF = 3
    DC_LDPC_BF = 7
    N_LDPC_BF_TARGET = 1000 
    # (Sua função LDPC_decode já lida com ajustes de N se necessário)

    print(f"INICIANDO SIMULAÇÃO LDPC ({DV_LDPC_BF},{DC_LDPC_BF}) N~{N_LDPC_BF_TARGET} (BitFlip-BSC)")
    start_time_ldpc_bf = time.time()
    try:
        ber_ldpc_bitflip = LDPC_decode( # Sua função para bit-flipping sobre BSC
            dc=DC_LDPC_BF,
            dv=DV_LDPC_BF,
            possiveis_P=channel_p_values,
            target_N=N_LDPC_BF_TARGET - (N_LDPC_BF_TARGET % DC_LDPC_BF)
        )
        if ber_ldpc_bitflip:
            all_ber_results_for_plot[f"LDPC ({DV_LDPC_BF},{DC_LDPC_BF}) N~{N_LDPC_BF_TARGET} (BitFlip-BSC)"] = ber_ldpc_bitflip
            print(f"BER LDPC (BitFlip-BSC): {[f'{b:.3e}' for b in ber_ldpc_bitflip]}")
        else:
            print("Simulação LDPC (BitFlip-BSC) não produziu resultados BER.")
    except Exception as e:
        print(f"Erro durante a simulação LDPC (BitFlip-BSC): {e}")
    end_time_ldpc_bf = time.time()
    print(f"Simulação LDPC (BitFlip-BSC) concluída em {end_time_ldpc_bf - start_time_ldpc_bf:.2f} segundos.")
    print("-" * 50)

    # --- 4. Simulação LDPC com LLR e Propagação de Crença usando AWGN ---
    # (Usando sua função `run_simulation` que chama `LDPC_decode_llr` ou `LDPC_decode_llr_structured`)
    DV_LDPC_LLR = 3
    DC_LDPC_LLR = 7
    N_LDPC_LLR_TARGET = 1000
    MAX_DECODER_ITER_LDPC_LLR = 50
    NUM_CODEWORDS_LDPC_LLR = 1000 # Ajuste para tempo de simulação vs precisão

    print(f"INICIANDO SIMULAÇÃO LDPC ({DV_LDPC_LLR},{DC_LDPC_LLR}) N~{N_LDPC_LLR_TARGET} (LLR-AWGN)")
    start_time_ldpc_llr = time.time()
    try:
        # Certifique-se que `run_simulation` e `LDPC_decode_llr` (ou _structured) estão definidas
        # e usam o canal AWGN (com `set_gaussian_channel_values`).
        ber_ldpc_llr, ebn0_axis_llr = run_simulation(
            dv=DV_LDPC_LLR, dc=DC_LDPC_LLR, N_target=N_LDPC_LLR_TARGET,
            eb_n0_db_range=EB_N0_DB_RANGE,
            max_iterations_decoder=MAX_DECODER_ITER_LDPC_LLR,
            num_codewords_to_simulate=NUM_CODEWORDS_LDPC_LLR,
            # min_errors_to_collect=100, # Adicione se quiser e se run_simulation suportar
            amplitude_bpsk=AMPLITUDE_BPSK
        )
        if ber_ldpc_llr:
            all_ber_results_for_plot[f"LDPC ({DV_LDPC_LLR},{DC_LDPC_LLR}) N~{N_LDPC_LLR_TARGET} (LLR-AWGN)"] = ber_ldpc_llr
            print(f"BER LDPC (LLR-AWGN): {[f'{b:.3e}' for b in ber_ldpc_llr]}")
        else:
            print("Simulação LDPC (LLR-AWGN) não produziu resultados BER.")
    except NameError as e:
        print(f"Erro: Função 'run_simulation' ou suas dependências (como LDPC_decode_llr) não definidas? {e}")
    except Exception as e:
        print(f"Erro durante a simulação LDPC (LLR-AWGN): {e}")
    end_time_ldpc_llr = time.time()
    print(f"Simulação LDPC (LLR-AWGN) concluída em {end_time_ldpc_llr - start_time_ldpc_llr:.2f} segundos.")
    print("-" * 50)

    # --- 5. Simulação LDPC (3,6) com LLR sobre AWGN ---
    DV_LDPC_36 = 3
    DC_LDPC_36 = 6
    N_LDPC_36 = 1000

    print(f"INICIANDO SIMULAÇÃO LDPC ({DV_LDPC_36},{DC_LDPC_36}) N~{N_LDPC_36} (LLR-AWGN)")
    start_time_ldpc_36 = time.time()
    try:
        ber_ldpc_36, ebn0_axis_ldpc_36 = run_simulation(
            dv=DV_LDPC_36,
            dc=DC_LDPC_36,
            N_target=N_LDPC_36,
            eb_n0_db_range=EB_N0_DB_RANGE,
            max_iterations_decoder=MAX_DECODER_ITER_LDPC_LLR,
            num_codewords_to_simulate=NUM_CODEWORDS_LDPC_LLR,
            amplitude_bpsk=AMPLITUDE_BPSK
        )
        if ber_ldpc_36:
            all_ber_results_for_plot[f"LDPC ({DV_LDPC_36},{DC_LDPC_36}) N~{N_LDPC_36} (LLR-AWGN)"] = ber_ldpc_36
            print(f"BER LDPC (LLR-AWGN, dv=3, dc=6): {[f'{b:.3e}' for b in ber_ldpc_36]}")
        else:
            print("Simulação LDPC (3,6) não produziu resultados BER.")
    except Exception as e:
        print(f"Erro durante a simulação LDPC (3,6): {e}")
    end_time_ldpc_36 = time.time()
    print(f"Simulação LDPC (3,6) concluída em {end_time_ldpc_36 - start_time_ldpc_36:.2f} segundos.")
    print("-" * 50)
    
    # --- Plotagem de Todos os Resultados ---
    plot_ber_comparison(
        eb_n0_db_axis=EB_N0_DB_RANGE, # Usamos este como eixo X comum
        ber_data_dict=all_ber_results_for_plot,
        filename=f"comparacao_TODOS_LDPC_Hamming.png"
    )

    print("Todas as simulações e plotagem concluídas.")
    
main()