import random
import matplotlib.pyplot as plt
from Lab_1.lab1 import hamming_simulation 
import os
import time
import csv
import math
import numpy as np
class VNode:
    def __init__(self):
        self.val_r = 0
        self.llr_channel = 0
        self.c_nodes = []
        self.mensagens_entrada_llr = [] 
        self.mensagens_saida_llr = []   
        self.dv = 0

    def add_c_node(self, cnode):
        if not isinstance(cnode, CNode):
            raise TypeError("Esperado um objeto do tipo CNode")
        self.c_nodes.append(cnode)
        self.dv = len(self.c_nodes)

    def initialize_llr_storage(self):
        """Chamado após todas as conexões serem feitas para alocar as listas de mensagens."""
        self.mensagens_entrada_llr = [0.0] * self.dv
        self.mensagens_saida_llr = [0.0] * self.dv
class CNode:
    def __init__(self):
        self.v_nodes = []
        self.mensagens_entrada_llr = [] 
        self.mensagens_saida_llr = []   
        self.dc = 0
        

    def add_v_node(self, vnode):
        if not isinstance(vnode, VNode):
            raise TypeError("Esperado um objeto do tipo VNode")
        self.v_nodes.append(vnode)
        self.dc = len(self.v_nodes)

    def initialize_llr_storage(self):
        """Chamado após todas as conexões serem feitas para alocar as listas de mensagens."""
        self.mensagens_entrada_llr = [0.0] * self.dc
        self.mensagens_saida_llr = [0.0] * self.dc

def LDPC(dv, dc, N):
    N = N - (N % dc)
    if (N * dv) % dc != 0:
        raise ValueError(f"N * dv ({N*dv}) deve ser divisível por dc ({dc}) para um grafo LDPC regular.")
    if N <= 0 or dv <= 0 or dc <= 0: 
        raise ValueError("N, dv e dc devem ser inteiros positivos.")

    M = (N * dv) // dc
    if M == 0:
        print("Aviso: Parâmetros resultam em M=0 check nodes. Retornando listas vazias.")
        return [], []


    all_vnodes = [VNode() for _ in range(N)]
    all_cnodes = [CNode() for _ in range(M)]

    # 1) Cria lista de stubs: cada v-node i (índice) aparece dv vezes
    stubs = []
    for v_idx in range(N):
        stubs.extend([v_idx] * dv) 

    random.shuffle(stubs)

    # 2) Para cada c-node, seleciona dc stubs distintos
    edges_to_place_total = N * dv
    edges_placed_count = 0

    for c_node_idx in range(M): 
        cnode_current = all_cnodes[c_node_idx]
        chosen_v_indices_for_this_cnode = set() 

        for _ in range(dc):
            if not stubs:
                raise RuntimeError(f"Faltaram stubs ao tentar conectar c-node {c_node_idx}. "
                                   f"Total de arestas a colocar: {edges_to_place_total}, colocadas: {edges_placed_count}. "
                                   "Isso indica um problema na lógica ou parâmetros extremos.")

            attempts_for_current_edge = 0       
            max_attempts_heuristic = max(len(stubs) * 2, N, 100)

            while True:
                attempts_for_current_edge += 1
                if attempts_for_current_edge > max_attempts_heuristic:
                    raise RuntimeError(
                        f"Excedido o máximo de {max_attempts_heuristic} tentativas para encontrar um v-node distinto "
                        f"para o c-node {c_node_idx}. Stubs restantes: {len(stubs)}. "
                        f"V-nodes já escolhidos para este c-node: {len(chosen_v_indices_for_this_cnode)}/{dc}. "
                        "Considere verificar os parâmetros (dv, dc, N) ou aumentar a heurística de max_attempts."
                    )

                if not stubs:
                     raise RuntimeError(f"Stubs esgotados inesperadamente dentro do loop de tentativas para c-node {c_node_idx}.")

                idx_in_stubs_list = random.randrange(len(stubs))
                v_node_idx_of_stub = stubs[idx_in_stubs_list] 

                if v_node_idx_of_stub not in chosen_v_indices_for_this_cnode:
                    chosen_v_indices_for_this_cnode.add(v_node_idx_of_stub)

                    vnode_to_connect = all_vnodes[v_node_idx_of_stub]

                    # Conecta
                    cnode_current.add_v_node(vnode_to_connect)
                    vnode_to_connect.add_c_node(cnode_current)
                    edges_placed_count += 1

                    # Remove o stub usado de forma eficiente (O(1))
                    # Trocando com o último elemento e depois pop
                    stubs[idx_in_stubs_list] = stubs[-1]
                    stubs.pop()
                    break 

    if stubs: 
        print(f"Aviso: {len(stubs)} stubs restantes após a construção do grafo. Deveria ser 0.")

    save_graph_to_csv(all_vnodes, all_cnodes)

    return all_vnodes, all_cnodes

def save_graph_to_csv(all_vnodes, all_cnodes, filename="ldpc_graph_versao2.csv"):
    """
    Salva a estrutura do grafo LDPC em um arquivo CSV com índices baseados em 1.
    O arquivo CSV é salvo no mesmo diretório que o script Python que contém esta função.
    A linha física `k` do arquivo (onde a primeira linha é k=1) corresponde ao v-node `all_vnodes[k-1]`.
    A linha conterá os ÍNDICES BASEADOS EM 1 dos c-nodes conectados,
    obtidos pela busca desses c-nodes na lista `all_cnodes` e adicionando 1.

    Args:
        all_vnodes (list): Lista de objetos VNode.
        all_cnodes (list): Lista de objetos CNode. Essencial para encontrar os índices.
        filename (str): Nome base do arquivo CSV a ser salvo (ex: "ldpc_graph.csv").
    """
    if not all_vnodes:
        print("Nenhum v-node para salvar.")
        return
    if not all_cnodes and any(v.c_nodes for v in all_vnodes): 
        print("Nenhum c-node fornecido, mas os v-nodes têm conexões. Não é possível determinar os índices dos c-nodes.")
        print("Não é possível salvar o grafo no formato desejado sem a lista de c-nodes.")
        return

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
        print(f"Aviso: __file__ não definido. Salvando '{filename}' no diretório de trabalho atual: {script_dir}")


    full_path = os.path.join(script_dir, filename)

    with open(full_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for vnode_obj in all_vnodes:
            connected_cnode_indices_one_based = []
            
            for cnode_connection in vnode_obj.c_nodes:
                try:
                    zero_based_idx = all_cnodes.index(cnode_connection)
                    one_based_idx = zero_based_idx + 1
                    connected_cnode_indices_one_based.append(one_based_idx)
                except ValueError:
                    print(f"Alerta: Objeto CNode conectado a um VNode não foi encontrado na lista principal 'all_cnodes'. "
                          f"Isso não deveria acontecer. Objeto CNode: {cnode_connection}")
                    
            writer.writerow(sorted(connected_cnode_indices_one_based))

    print(f"Grafo salvo em {full_path} com índices baseados em 1.")

def set_gaussian_channel_values(all_vnodes, transmitted_codeword_bits, N0, amplitude=1.0):
    """
    Simula o canal AWGN para uma dada palavra código transmitida e calcula os LLRs do canal.

    Args:
        all_vnodes (list): Uma lista de objetos VNode.
        transmitted_codeword_bits (list of int): Lista de bits (0 ou 1) da palavra código transmitida.
        N0 (float): Densidade espectral de potência do ruído (unilateral).
        amplitude (float): Amplitude do sinal BPSK (e.g., 1.0).
    """
    if len(all_vnodes) != len(transmitted_codeword_bits):
        raise ValueError("Número de vnodes deve ser igual ao tamanho da palavra código transmitida.")

    variance = N0 / 2.0
    #Isso é realmente necessário?
    if variance <= 0: 
        print("Aviso: Variância do canal <= 0. LLRs podem ser instáveis.")
        for i, vnode in enumerate(all_vnodes):
            transmitted_symbol = amplitude if transmitted_codeword_bits[i] == 0 else -amplitude
            vnode.val_r = transmitted_symbol
            vnode.llr_channel = 100.0 * np.sign(transmitted_symbol) if transmitted_symbol != 0 else 0.0
        return

    std_dev = math.sqrt(variance)

    for i, vnode in enumerate(all_vnodes):
        mean_signal = amplitude if transmitted_codeword_bits[i] == 0 else -amplitude
        
        vnode.val_r = random.gauss(mean_signal, std_dev) 
        vnode.llr_channel = (2 * vnode.val_r * amplitude) / variance

def LDPC_decode_llr(all_vnodes, all_cnodes, initial_channel_llrs, max_iter):
    """
    Decodifica uma palavra código LDPC seguindo estritamente a sequência do "Algoritmo",
    usando o critério de parada "produto dos sinais" nos C-nodes.

    Args:
        all_vnodes (list): Lista de objetos VNode do grafo.
        all_cnodes (list): Lista de objetos CNode do grafo.
        initial_channel_llrs (list): LLRs do canal para cada v-node.
        max_iter (int): Número máximo de iterações.

    Returns:
        list of int: A palavra código decodificada (0s e 1s).
    """
    N = len(all_vnodes)
    if N == 0: return []
    if len(initial_channel_llrs) != N:
        raise ValueError("Número de LLRs do canal não corresponde ao número de V-nodes.")

    for i, vnode in enumerate(all_vnodes):
        vnode.llr_channel = initial_channel_llrs[i]

    # Passo 3 do Algoritmo: Inicie o contador de iterações.
    for iteration in range(max_iter):
        # --- Passo 4 do Algoritmo: Cálculo das mensagens de SAÍDA dos V-nodes (L_{v->c}) ---
        # L_{v->c_j} = Lc(v) + sum(L_{c_k->v} para k != j)
        for vnode in all_vnodes:
            if vnode.dv == 0: continue
            for c_target_idx in range(vnode.dv): 
                sum_llr_from_other_cnodes = 0.0
                for c_source_idx in range(vnode.dv): 
                    if c_source_idx == c_target_idx: 
                        continue
                    sum_llr_from_other_cnodes += vnode.mensagens_entrada_llr[c_source_idx]
                vnode.mensagens_saida_llr[c_target_idx] = vnode.llr_channel + sum_llr_from_other_cnodes

        # --- Passo 5 do Algoritmo: Teste se o critério de parada foi atingido ---
        # Critério 2: "para todos c_nodes, o produto de todas as mensagens chegando no node tem que ser positivo."    
        all_cnodes_ok_by_product_of_signs = True
        if not all_cnodes: 
            all_cnodes_ok_by_product_of_signs = False 
        
        for cnode_check in all_cnodes:
            if cnode_check.dc == 0: 
                continue
            
            product_of_incoming_signs = 1.0
            for v_node_source in cnode_check.v_nodes:
                try:
                    idx_of_cnode_in_vsource = v_node_source.c_nodes.index(cnode_check)
                    llr_v_to_c = v_node_source.mensagens_saida_llr[idx_of_cnode_in_vsource]
                    product_of_incoming_signs *= np.sign(llr_v_to_c) if llr_v_to_c != 0 else 1.0
                except ValueError:
                    raise RuntimeError(f"Consistência do grafo: CNode {cnode_check} não encontrado em VNode {v_node_source}")
            
            if product_of_incoming_signs < 0: 
                all_cnodes_ok_by_product_of_signs = False
                break 

        if all_cnodes_ok_by_product_of_signs and all_cnodes: 
            decoded_bits = [0] * N
            for i_vn, vn_decision in enumerate(all_vnodes):
                sum_incoming_c_llrs_for_decision = sum(vn_decision.mensagens_entrada_llr) 
                app_llr = vn_decision.llr_channel + sum_incoming_c_llrs_for_decision
                decoded_bits[i_vn] = 0 if app_llr >= 0 else 1
            
            is_codeword_truly_valid = True
            if not all_cnodes: is_codeword_truly_valid = False 

            for cnode_verify in all_cnodes:
                current_parity_sum = 0
                for v_conn_verify in cnode_verify.v_nodes:
                    idx_global_v = all_vnodes.index(v_conn_verify)
                    current_parity_sum = (current_parity_sum + decoded_bits[idx_global_v]) % 2
                if current_parity_sum != 0:
                    is_codeword_truly_valid = False
                    break
            
            if is_codeword_truly_valid:
                return decoded_bits

        # Se o critério de parada NÃO foi atingido (ou a palavra decodificada não era válida):
        # --- Transferência das L_{v->c} para as entradas dos C-nodes ---
        # (Necessário para o Passo 6)
        # As vnode.mensagens_saida_llr (L_{v->c} do Passo 4) tornam-se cnode.mensagens_entrada_llr.
        for vnode_transfer in all_vnodes:
            for c_conn_idx, target_cnode in enumerate(vnode_transfer.c_nodes):
                try:
                    v_idx_in_cnode = target_cnode.v_nodes.index(vnode_transfer)
                    target_cnode.mensagens_entrada_llr[v_idx_in_cnode] = vnode_transfer.mensagens_saida_llr[c_conn_idx]
                except ValueError:
                     raise RuntimeError(f"Consistência do grafo: VNode {vnode_transfer} não encontrado em CNode {target_cnode}")

        # --- Passo 6 do Algoritmo: Cálculo das mensagens de SAÍDA dos C-nodes (L_{c->v}) ---
        # L_{c->v_i} ≈ min_{k!=i} |L_{v_k->c}| * product_{k!=i} sign(L_{v_k->c})
        for cnode in all_cnodes:
            if cnode.dc == 0: continue
            for v_target_idx in range(cnode.dc):
                min_abs_llr_others = float('inf')
                prod_signs_others = 1.0
                for v_source_idx in range(cnode.dc): 
                    if v_source_idx == v_target_idx: 
                        continue
                    incoming_llr = cnode.mensagens_entrada_llr[v_source_idx]
                    min_abs_llr_others = min(min_abs_llr_others, abs(incoming_llr))
                    prod_signs_others *= np.sign(incoming_llr) if incoming_llr != 0 else 1.0
                
                if cnode.dc == 1: 
                    cnode.mensagens_saida_llr[v_target_idx] = 0.0
                elif cnode.dc == 2: 
                    idx_other = 1 - v_target_idx
                    cnode.mensagens_saida_llr[v_target_idx] = cnode.mensagens_entrada_llr[idx_other]
                else: 
                    cnode.mensagens_saida_llr[v_target_idx] = prod_signs_others * min_abs_llr_others
        
        # --- Transferência das L_{c->v} para as entradas dos V-nodes ---
        # (Necessário para o Passo 4 da PRÓXIMA iteração)
        # As cnode.mensagens_saida_llr (L_{c->v} do Passo 6) tornam-se vnode.mensagens_entrada_llr.
        for cnode_transfer in all_cnodes:
            for v_conn_idx, target_vnode in enumerate(cnode_transfer.v_nodes):
                try:
                    c_idx_in_vnode = target_vnode.c_nodes.index(cnode_transfer)
                    target_vnode.mensagens_entrada_llr[c_idx_in_vnode] = cnode_transfer.mensagens_saida_llr[v_conn_idx]
                except ValueError:
                    raise RuntimeError(f"Consistência do grafo: CNode {cnode_transfer} não encontrado em VNode {target_vnode}")

        # Passo 7 do Algoritmo: "Incremente ... e retorne ao passo 4" é feito pelo loop `for`.

    # --- Passo 8 do Algoritmo: Decida sobre os bits transmitidos (após max_iter) ---
    final_decoded_bits = [0] * N
    for i_vn_final, vn_final_decision in enumerate(all_vnodes):
        sum_incoming_c_llrs_final = sum(vn_final_decision.mensagens_entrada_llr)
        app_llr_final = vn_final_decision.llr_channel + sum_incoming_c_llrs_final
        final_decoded_bits[i_vn_final] = 0 if app_llr_final >= 0 else 1
        
    return final_decoded_bits

def eb_n0_dB_to_N0(eb_n0_db_value, code_rate, amplitude=1.0):
    """
    Converte Eb/N0 em dB para o valor linear de N0.
    Assume BPSK onde Eb = A^2 * Tb. Se Tb=1, Eb = A^2.
    Para códigos, usamos Eb_channel_bit / N0.
    Se o code_rate R = k/n, então Eb_info_bit = Eb_channel_bit / R.
    A simulação geralmente trabalha com Eb_channel_bit / N0.
    Se a entrada eb_n0_db_value é para bits de informação, ajuste é necessário.
    Assumindo aqui que eb_n0_db_value é para bits de canal (codificados).

    Args:
        eb_n0_db_value (float): Valor de Eb/N0 em dB.
        code_rate (float): Taxa do código (k/N). Para LDPC (dv,dc), R = 1 - dv/dc.
        amplitude (float): Amplitude do sinal BPSK.

    Returns:
        float: Valor de N0 (densidade espectral de potência do ruído unilateral).
    """
    eb_n0_linear = 10**(eb_n0_db_value / 10.0)
    eb_channel_n0_linear = eb_n0_linear * code_rate 
    
    if eb_channel_n0_linear == 0: 
        return float('inf') 
    
    N0 = amplitude**2 / eb_channel_n0_linear
    return N0

def run_simulation(dv, dc, N_target, eb_n0_db_range, max_iterations_decoder, num_codewords_to_simulate, amplitude_bpsk=1.0): # Removido min_errors_to_collect
    """
    Executa a simulação de Monte Carlo para um código LDPC.
    Simula um número fixo de palavras-código por ponto Eb/N0.
    """
    
    print(f"### Iniciando Simulação de Monte Carlo para LDPC ({dv},{dc}), N_alvo={N_target} ###")
    print(f"### Simulando {num_codewords_to_simulate} palavras-código por ponto Eb/N0 ###") # Mensagem atualizada
    
    N_actual = N_target - (N_target % dc)
    if N_actual == 0:
        print("N_actual resultou em 0. Simulação abortada.")
        return [], []
        
    if dc == 0:
        print("Erro: dc (grau do C-node) não pode ser zero. Simulação abortada.")
        return [],[]
    
    code_rate_estimate = 1.0 - (dv / dc)
    if not (0 < code_rate_estimate < 1):
        M_actual = (N_actual * dv) // dc
        K_actual = N_actual - M_actual
        if N_actual > 0 and K_actual > 0:
            code_rate_estimate = K_actual / N_actual
            print(f"Aviso: dv={dv}, dc={dc} resulta em taxa fora do comum. Taxa K/N estimada: {code_rate_estimate:.3f}")
        else:
            print(f"Alerta Crítico: dv={dv}, dc={dc}, N={N_actual} resulta em K<=0. Usando taxa nominal de 0.5.")
            code_rate_estimate = 0.5 
    print(f"Taxa do código (R) estimada: {code_rate_estimate:.3f}")

    print(f"Construindo grafo LDPC para N={N_actual} (dv={dv}, dc={dc})...")
    start_graph_time = time.time()
    all_vnodes, all_cnodes = LDPC(dv, dc, N_actual)
    end_graph_time = time.time()
    if not all_vnodes:
        print(f"Falha ao gerar o grafo LDPC para N={N_actual}. Simulação abortada.")
        return [], []
    print(f"Grafo construído em {end_graph_time - start_graph_time:.2f}s. Inicializando LLR storage...")
    
    for vnode in all_vnodes: vnode.initialize_llr_storage()
    for cnode in all_cnodes: cnode.initialize_llr_storage()
    print("LLR storage inicializado.")

    ber_results = []
    transmitted_codeword = [0] * N_actual

    print("\n--- Iniciando Loop de Eb/N0 ---")
    for point_idx, eb_n0_db in enumerate(eb_n0_db_range):
        N0 = eb_n0_dB_to_N0(eb_n0_db, code_rate_estimate, amplitude_bpsk)
        print(f"\n[Ponto Eb/N0 {point_idx+1}/{len(eb_n0_db_range)}] Eb/N0 = {eb_n0_db:.2f} dB (N0 = {N0:.3e}, Variância do ruído = {N0/2:.3e})")

        total_bits_simulated_this_point = 0
        total_bit_errors_this_point = 0
        
        start_time_ebn0_point = time.time()
        for i_codeword in range(num_codewords_to_simulate): # Loop fixo
            set_gaussian_channel_values(all_vnodes, transmitted_codeword, N0, amplitude_bpsk)
            initial_llrs_for_decoder = [v.llr_channel for v in all_vnodes]
            decoded_codeword = LDPC_decode_llr(all_vnodes, all_cnodes, initial_llrs_for_decoder, max_iterations_decoder)

            current_errors = 0
            for bit_idx in range(N_actual):
                if transmitted_codeword[bit_idx] != decoded_codeword[bit_idx]:
                    current_errors += 1
            
            total_bit_errors_this_point += current_errors
            total_bits_simulated_this_point += N_actual

            if (i_codeword + 1) % (max(1, num_codewords_to_simulate // 20)) == 0 or (i_codeword + 1) == num_codewords_to_simulate:
                 current_ber_estimate = total_bit_errors_this_point / total_bits_simulated_this_point if total_bits_simulated_this_point > 0 else 0
                 elapsed_time_this_ebn0 = time.time() - start_time_ebn0_point
                 print(f"    Palavra {i_codeword+1:6d}/{num_codewords_to_simulate:6d} | Erros Bit: {total_bit_errors_this_point:6d} | Total Bits: {total_bits_simulated_this_point:9d} | BER Est.: {current_ber_estimate:.2e} | Tempo: {elapsed_time_this_ebn0:.1f}s")

        
        end_time_ebn0_point = time.time()
        ber = total_bit_errors_this_point / total_bits_simulated_this_point if total_bits_simulated_this_point > 0 else 0.0
        ber_results.append(ber)
        print(f"  Resultado Final para Eb/N0 = {eb_n0_db:.1f} dB: BER = {ber:.4e} (Total Erros: {total_bit_errors_this_point}, Total Bits: {total_bits_simulated_this_point})")
        print(f"  Tempo para este ponto Eb/N0: {end_time_ebn0_point - start_time_ebn0_point:.2f} segundos.")

    print("\n### Simulação de Monte Carlo Concluída ###")
    return ber_results, eb_n0_db_range

def plot_ber_vs_ebn0(eb_n0_db_values, ber_values, code_label, filename="ber_ldpc_plot.png"):
    """
    Plota BER vs Eb/N0.
    """
    plt.figure(figsize=(10, 7))
    plt.semilogy(eb_n0_db_values, ber_values, marker='o', linestyle='-', label=code_label)
    
    eb_n0_linear_uncoded = [10**(x/10.0) for x in eb_n0_db_values]
    ber_uncoded_bpsk = [0.5 * math.erfc(math.sqrt(x)) for x in eb_n0_linear_uncoded] # Q(x) = 0.5 * erfc(x/sqrt(2))
    plt.semilogy(eb_n0_db_values, ber_uncoded_bpsk, marker='s', linestyle='--', label="BPSK Não Codificado (Teórico)")

    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("Bit Error Rate (BER)")
    plt.title(f"Desempenho LDPC: {code_label}")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.ylim(1e-6, 1e0) 
    plt.xlim(min(eb_n0_db_values), max(eb_n0_db_values))

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
    full_path = os.path.join(script_dir, filename)
    plt.savefig(full_path)
    print(f"Gráfico salvo em: {full_path}")
    plt.show() 


# --- Configurações da Simulação ---
if __name__ == "__main__":
    DV = 3
    DC = 7
    N = 1000
    EB_N0_DB_RANGE = np.arange(0.0, 5.1, 0.5).tolist()
    MAX_DECODER_ITERATIONS = 50
    NUM_CODEWORDS_PER_POINT = 10000 
    AMPLITUDE = 1.0

    print("-" * 50)
    print(f"INICIANDO SIMULAÇÃO PARA N_TARGET = {N}")
    start_time = time.time()
    ber, ebn0 = run_simulation(
        dv=DV, dc=DC, N_target=N,
        eb_n0_db_range=EB_N0_DB_RANGE,
        max_iterations_decoder=MAX_DECODER_ITERATIONS,
        num_codewords_to_simulate=NUM_CODEWORDS_PER_POINT,
        amplitude_bpsk=AMPLITUDE
    )
    end_time = time.time()
    print(f"Simulação para N_target={N} concluída em {end_time - start_time:.2f} segundos.")
    
    if ber: 
        plot_ber_vs_ebn0(
            ebn0, ber,
            code_label=f"LDPC ({DV},{DC}), N~{N}",
            filename=f"ber_ldpc_N{N}_dv{DV}_dc{DC}_versao2.png"
        )

    print("-" * 50)
    print("Todas as simulações concluídas.")
    

