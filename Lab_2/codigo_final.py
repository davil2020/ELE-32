import random
import matplotlib.pyplot as plt
from Lab_1.lab1 import hamming_simulation 
import os
import time
import csv


class VNode:
    def __init__(self):
        self.val = 0
        self.c_nodes = []
        self.dv = 0
        self.count_odd_c_node = 0

    def add_c_node(self, cnode):
        if not isinstance(cnode, CNode):
            raise TypeError("Esperado um objeto do tipo CNode")
        self.c_nodes.append(cnode)
        self.dv = len(self.c_nodes)


class CNode:
    def __init__(self):
        self.v_nodes = []
        self.dc = 0
        self.parid_check = 0

    def add_v_node(self, vnode):
        if not isinstance(vnode, VNode):
            raise TypeError("Esperado um objeto do tipo VNode")
        self.v_nodes.append(vnode)
        self.dc = len(self.v_nodes)

def LDPC(dv, dc, N):
    if (N * dv) % dc != 0:
        # Este erro é fundamental, a geração não pode prosseguir
        raise ValueError(f"N * dv ({N*dv}) deve ser divisível por dc ({dc}) para um grafo LDPC regular.")
    if N <= 0 or dv <= 0 or dc <= 0: # Checa se são positivos
        raise ValueError("N, dv e dc devem ser inteiros positivos.")

    M = (N * dv) // dc
    if M == 0: # Caso onde N*dv < dc, o que significa nenhum c-node
        print("Aviso: Parâmetros resultam em M=0 check nodes. Retornando listas vazias.")
        return [], []


    all_vnodes = [VNode() for _ in range(N)]
    all_cnodes = [CNode() for _ in range(M)]

    # 1) Cria lista de stubs: cada v-node i (índice) aparece dv vezes
    stubs = []
    for v_idx in range(N):
        stubs.extend([v_idx] * dv) # .extend é um pouco mais eficiente que laço com += lista

    # Opcional: Embaralhar os stubs uma vez no início.
    # Pode ajudar a quebrar padrões se a ordem inicial dos stubs tiver algum viés,
    # mas com random.randrange, o impacto na performance/aleatoriedade final é mínimo.
    random.shuffle(stubs)

    # 2) Para cada c-node, seleciona dc stubs distintos
    edges_to_place_total = N * dv
    edges_placed_count = 0

    for c_node_idx in range(M): # Itera pelos índices dos c-nodes
        cnode_current = all_cnodes[c_node_idx]
        chosen_v_indices_for_this_cnode = set() # V-nodes já conectados a ESTE c-node

        # Cada c-node precisa de 'dc' conexões distintas
        for _ in range(dc):
            if not stubs:
                # Isso não deveria acontecer se N*dv == M*dc e tudo correu bem
                # Indica que ficamos sem stubs antes de completar as conexões.
                raise RuntimeError(f"Faltaram stubs ao tentar conectar c-node {c_node_idx}. "
                                   f"Total de arestas a colocar: {edges_to_place_total}, colocadas: {edges_placed_count}. "
                                   "Isso indica um problema na lógica ou parâmetros extremos.")

            attempts_for_current_edge = 0
            # A heurística para max_attempts pode ser importante para evitar loops infinitos
            # em casos patológicos (raros com parâmetros LDPC típicos).
            # Se len(stubs) fica muito pequeno e todos os stubs restantes são para v_nodes já em
            # chosen_v_indices_for_this_cnode, o loop pode demorar.
            max_attempts_heuristic = max(len(stubs) * 2, N, 100) # Ajustado para ser mais robusto

            while True:
                attempts_for_current_edge += 1
                if attempts_for_current_edge > max_attempts_heuristic:
                    raise RuntimeError(
                        f"Excedido o máximo de {max_attempts_heuristic} tentativas para encontrar um v-node distinto "
                        f"para o c-node {c_node_idx}. Stubs restantes: {len(stubs)}. "
                        f"V-nodes já escolhidos para este c-node: {len(chosen_v_indices_for_this_cnode)}/{dc}. "
                        "Considere verificar os parâmetros (dv, dc, N) ou aumentar a heurística de max_attempts."
                    )

                if not stubs: # Checagem extra dentro do loop de tentativas
                     raise RuntimeError(f"Stubs esgotados inesperadamente dentro do loop de tentativas para c-node {c_node_idx}.")

                # Escolhe um índice aleatório da lista de stubs
                idx_in_stubs_list = random.randrange(len(stubs))
                v_node_idx_of_stub = stubs[idx_in_stubs_list] # Este é o índice do v-node

                # Garante que este v-node ainda não foi conectado a ESTE c-node
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
                    break # Passa para a próxima conexão deste c-node
                # Se v_node_idx_of_stub JÁ ESTÁ em chosen_v_indices_for_this_cnode,
                # o loop while continua para tentar outro stub.

    if stubs: # Ao final, a lista de stubs deveria estar vazia
        # Isso é um erro de lógica se N*dv == M*dc
        print(f"Aviso: {len(stubs)} stubs restantes após a construção do grafo. Deveria ser 0.")
        # Poderia levantar um erro aqui também, pois indica um descompasso.
        # raise RuntimeError(f"{len(stubs)} stubs restantes. Total de arestas esperado: {edges_to_place_total}, colocadas: {edges_placed_count}")


    # Verificação final de graus (opcional, mas bom para depuração)
    # for i, vnode in enumerate(all_vnodes):
    #     if vnode.dv != dv:
    #         print(f"Alerta: VNode {i} tem grau {vnode.dv}, esperado {dv}")
    # for i, cnode in enumerate(all_cnodes):
    #     if cnode.dc != dc:
    #         print(f"Alerta: CNode {i} tem grau {cnode.dc}, esperado {dc}")


    return all_vnodes, all_cnodes

def save_graph_to_csv(all_vnodes, all_cnodes, filename="ldpc_graph.csv"):
    """
    Salva a estrutura do grafo LDPC em um arquivo CSV com índices baseados em 1.
    A linha física `k` do arquivo (onde a primeira linha é k=1) corresponde ao v-node `all_vnodes[k-1]`.
    A linha conterá os ÍNDICES BASEADOS EM 1 dos c-nodes conectados,
    obtidos pela busca desses c-nodes na lista `all_cnodes` e adicionando 1.

    Args:
        all_vnodes (list): Lista de objetos VNode.
        all_cnodes (list): Lista de objetos CNode. Essencial para encontrar os índices.
        filename (str): Nome do arquivo CSV a ser salvo.
    """
    if not all_vnodes:
        print("Nenhum v-node para salvar.")
        return
    if not all_cnodes and any(v.c_nodes for v in all_vnodes): # Se há v-nodes com conexões mas não c-nodes
        print("Nenhum c-node fornecido, mas os v-nodes têm conexões. Não é possível determinar os índices dos c-nodes.")
        print("Não é possível salvar o grafo no formato desejado sem a lista de c-nodes.")
        return

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Para cada objeto vnode na lista all_vnodes...
        # A k-ésima iteração aqui (começando em k=0) corresponderá à (k+1)-ésima linha no arquivo.
        for vnode_obj in all_vnodes:
            connected_cnode_indices_one_based = []
            # Para cada objeto cnode ao qual este vnode_obj está conectado...
            for cnode_connection in vnode_obj.c_nodes:
                try:
                    # Encontra o índice baseado em 0 do objeto cnode_connection na lista all_cnodes.
                    zero_based_idx = all_cnodes.index(cnode_connection)
                    # Adiciona 1 para torná-lo baseado em 1
                    one_based_idx = zero_based_idx + 1
                    connected_cnode_indices_one_based.append(one_based_idx)
                except ValueError:
                    print(f"Alerta: Objeto CNode conectado a um VNode não foi encontrado na lista principal 'all_cnodes'. "
                          f"Isso não deveria acontecer. CNode object: {cnode_connection}")
                    # Decide como lidar: pular, adicionar placeholder, etc.

            # Ordena os índices baseados em 1 para uma saída consistente
            writer.writerow(sorted(connected_cnode_indices_one_based))

    print(f"Grafo salvo em {filename} com índices baseados em 1.")

def BSC_bit_flip(all_vnodes, p):
    for vnode in all_vnodes:
        if random.random() < p:
            vnode.val ^= 1  # Flip: 0 vira 1 e 1 vira 0

def Lista_N(dc = 7):
    valores = [100, 200, 500, 1000]
    maiores_multiplos = [] 
    for v in valores:
        maior_multiplo = v - (v % dc)
        maiores_multiplos.append(maior_multiplo)
    return maiores_multiplos

def generate_vector_p():
    p = [0.5, 0.2, 0.1]
    q = [0.5, 0.2, 0.1]
    for i in range(1,5):
      for j in range(3):
        q[j] = q[j]/10
        p.append(q[j])
    return p

def cnodes_parid_check(all_cnodes):
    for cnode in all_cnodes:
        cnode.parid_check = 0
        for vnode in cnode.v_nodes:
            cnode.parid_check = (cnode.parid_check + vnode.val) % 2

def vnodes_bit_ajust(all_vnodes):
    is_ok = True
    vnode_to_ajust = all_vnodes[0]
    for vnode in all_vnodes:
        vnode.count_odd_c_node = 0
        for cnode in vnode.c_nodes:
            vnode.count_odd_c_node += cnode.parid_check
        if vnode.count_odd_c_node > 0:
            is_ok = False
        if vnode.count_odd_c_node > vnode_to_ajust.count_odd_c_node:
            vnode_to_ajust = vnode
    if not is_ok:
        vnode_to_ajust.val = (vnode_to_ajust.val + 1) % 2
    return is_ok

def LDPC_decode(dc, dv):
    # dc = 7
    # dv = 3
    possiveis_P = generate_vector_p()
    possiveis_N = Lista_N(dc)

    templates = {}
    for N in possiveis_N:
        print(f'templare relativo a {N} criado')
        templates[N] = LDPC(dv, dc, N)

    num_N = len(possiveis_N)
    num_P = len(possiveis_P)
    monte_carlo_runs = 1000

    # Inicializa os vetores de probabilidade para cada N
    prob_n0 = [0.0] * num_P
    prob_n1 = [0.0] * num_P
    prob_n2 = [0.0] * num_P
    prob_n3 = [0.0] * num_P

    prob_n = [prob_n0, prob_n1, prob_n2, prob_n3]

   
    

    for idx_p, p in enumerate(possiveis_P):
        for idx_n, N in enumerate(possiveis_N):
            total_bits = 0
            total_bits_errados = 0
            all_vnodes, all_cnodes = templates[N]

            for _ in range(monte_carlo_runs):
                for vnode in all_vnodes:
                    vnode.val = 0
                    vnode.count_odd_c_node = 0

                for cnode in all_cnodes:
                    cnode.parid_check = 0
                BSC_bit_flip(all_vnodes, p)
                max_iter = 50
                iteration = 0
                while iteration < max_iter:
                    # print([vnode.val for vnode in all_vnodes])
                    # print([cnode.parid_check for cnode in all_cnodes])
                    cnodes_parid_check(all_cnodes)
                    is_ok = vnodes_bit_ajust(all_vnodes)
                    if is_ok:
                        break
                    iteration += 1
                bits_errados = sum(vnode.val for vnode in all_vnodes)
                total_bits_errados += bits_errados
                total_bits += len(all_vnodes)
                # print([vnode.val for vnode in all_vnodes])
                # print([cnode.parid_check for cnode in all_cnodes])
                # print(iteration)

                
            # Calcula BER para este (p, N)
            prob_n[idx_n][idx_p] = total_bits_errados / total_bits

    total_N = sum(possiveis_N)
    prob_geral = []
    for idx_p in range(num_P):
        weighted_sum = (possiveis_N[0] * prob_n0[idx_p] +
                        possiveis_N[1] * prob_n1[idx_p] +
                        possiveis_N[2] * prob_n2[idx_p] +
                        possiveis_N[3] * prob_n3[idx_p])
        prob_geral.append(weighted_sum / total_N)
   
    
    return prob_n0, prob_n1, prob_n2, prob_n3, prob_geral

def plot_ldpc_results(prob_n0, prob_n1, prob_n2, prob_n3, prob_geral, dc, dv, possiveis_N, filename, error_hamming=None):
    vector_p = generate_vector_p()

    plt.figure(figsize=(10, 6))

    plt.semilogx(vector_p, prob_n0, marker='o', linestyle='-', label=f'LDPC N={possiveis_N[0]}')
    plt.semilogx(vector_p, prob_n1, marker='s', linestyle='--', label=f'LDPC N={possiveis_N[1]}')
    plt.semilogx(vector_p, prob_n2, marker='^', linestyle='-.', label=f'LDPC N={possiveis_N[2]}')
    plt.semilogx(vector_p, prob_n3, marker='v', linestyle=':', label=f'LDPC N={possiveis_N[3]}')
    plt.semilogx(vector_p, prob_geral, marker='x', linestyle='-', label='Média Ponderada')
    plt.semilogx(vector_p, vector_p, marker='d', linestyle=':', label='Não codificado')
    if error_hamming is not None:
        plt.semilogx(vector_p, error_hamming, marker='*', linestyle='--', label='Hamming (7,4)')
    plt.yscale('log')
    plt.gca().invert_xaxis()

    plt.xlabel("Probabilidade de erro (p)")
    plt.ylabel("Taxa de erro de bits")
    plt.title(f"Correção de Erros - LDPC (dv={dv}, dc={dc})")
    plt.grid(True, which="both", linestyle="--")
    plt.legend()

    # Salva o gráfico
    save_path = os.path.join(os.path.dirname(__file__), filename)
    plt.savefig(save_path)
    plt.close()


def main():
    configs = [
        (7, 3, "ldpc_dv3_dc7.png"),   # dv=3, dc=7
        (2, 1, "ldpc_dv1_dc2.png"),   # dv=1, dc=2
        (3, 2, "ldpc_dv2_dc3.png"),   # dv=2, dc=3
    ]

    hamming_probs = hamming_simulation()
    for dc, dv, filename in configs:
        possiveis_N = Lista_N(dc)
        print(f"\nRodando simulação para dv={dv}, dc={dc}")
        t0 = time.perf_counter()
        prob_n0, prob_n1, prob_n2, prob_n3, prob_geral = LDPC_decode(dc=dc, dv=dv)
        dt = time.perf_counter() - t0
        print(f"Tempo total: {dt:.2f}s")
        plot_ldpc_results(prob_n0, prob_n1, prob_n2, prob_n3, prob_geral, dc, dv, possiveis_N, filename, hamming_probs)
        print(f"Gráfico salvo como {filename}")

main()
# maior_multiplo = 1000 - (1000 % 7)
# [all_vnodes, all_cnodes]=LDPC(3,7,maior_multiplo)
# save_graph_to_csv(all_vnodes, all_cnodes)
# LDPC(3,7,98)
# print(Lista_N(dc=2))
# dc = 7 e dv=3

# [all_vnodes, all_cnodes]= LDPC(3,7,14)
# for i, vnode in enumerate(all_vnodes):
#     c_indices = [ all_cnodes.index(cnode) for cnode in vnode.c_nodes ]
#     print(f"VNode {i}: conectado a CNodes {c_indices}")
# print(all_vnodes)
# print(all_cnodes)


