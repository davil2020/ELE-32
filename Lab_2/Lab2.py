import random
import matplotlib.pyplot as plt
from Lab_1.lab1 import hamming_simulation 
import os
import time

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
        raise ValueError("N * dv deve ser divisível por dc para um LDPC regular.")

    M = (N * dv) // dc
    all_vnodes = [VNode() for _ in range(N)]
    all_cnodes = [CNode() for _ in range(M)]

    # 1) Cria lista de stubs: cada v-node i aparece dv vezes
    stubs = []
    for i in range(N):
        stubs += [i] * dv

    # 2) Para cada c-node, seleciona dc stubs distintos
    for cnode in all_cnodes:
        escolhidos = set()
        for _ in range(dc):
            while True:
                # escolhe stub aleatório
                idx = random.randrange(len(stubs))
                v_idx = stubs[idx]
                # garante distinctidade dentro deste c-node
                if v_idx not in escolhidos:
                    escolhidos.add(v_idx)
                    # conecta
                    vnode = all_vnodes[v_idx]
                    cnode.add_v_node(vnode)
                    vnode.add_c_node(cnode)
                    # remove stub em O(1)
                    stubs[idx] = stubs[-1]
                    stubs.pop()
                    break
                # caso contrário, repete até achar stub válido

    return all_vnodes, all_cnodes

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

def LDPC_decode(dc, dv, possiveis_P=None, target_N=None):
    # dc = 7
    # dv = 3
    print(f"DEBUG LDPC_decode Lab2: `possiveis_P` recebido: {possiveis_P} (Tipo: {type(possiveis_P)})") # D
    if possiveis_P is None:
        print("Nenhum `possiveis_P` fornecido, gerando um padrão...")
        possiveis_P = generate_vector_p()

    if target_N is not None:
        print(f"Simulando LDPC para N específico = {target_N}")
        N_sim = target_N 
        
        print(f"Construindo template para N={N_sim}...")
        all_vnodes, all_cnodes = LDPC(dv, dc, N_sim)
        
        if not all_vnodes:
            print(f"Não foi possível criar um grafo para N={N_sim} com dv={dv}, dc={dc}. Retornando lista de BER vazia.")
            return [] # Ou levante um erro, ou retorne None

        N_actual = len(all_vnodes)
        if N_actual == 0: # Checagem extra
            print(f"Template para N={N_sim} resultou em 0 vnodes. Retornando lista de BER vazia.")
            return []

        print(f"Template criado para N_actual={N_actual} (target_N={target_N}).")

        prob_para_target_N = [0.0] * len(possiveis_P)
        monte_carlo_runs = 1000 # Ou configure como parâmetro da função

        for idx_p, p_val in enumerate(possiveis_P):
            # print(f"  Simulando N={N_actual}, p={p_val:.5f}...")
            total_bits = 0
            total_bits_errados = 0

            for _ in range(monte_carlo_runs):
                # Resetar nós para cada run do Monte Carlo
                for vnode in all_vnodes: vnode.val = 0; vnode.count_odd_c_node = 0
                for cnode in all_cnodes: cnode.parid_check = 0
                
                BSC_bit_flip(all_vnodes, p_val)
                max_iter_decoder = 50 # Ou configure
                iteration = 0
                while iteration < max_iter_decoder:
                    cnodes_parid_check(all_cnodes)
                    is_ok = vnodes_bit_ajust(all_vnodes)
                    if is_ok: break
                    iteration += 1
                
                bits_errados_run = sum(vnode.val for vnode in all_vnodes)
                total_bits_errados += bits_errados_run
                total_bits += N_actual # Usar N_actual (len(all_vnodes))
            
            prob_para_target_N[idx_p] = total_bits_errados / total_bits if total_bits > 0 else 0.0
            # print(f"    BER para N={N_actual}, p={p_val:.5f}: {prob_para_target_N[idx_p]:.5e}")
        
        return prob_para_target_N
    else:
    
        possiveis_N = Lista_N(dc)

        templates = {}
        for N in possiveis_N:
            print(f'templare relativo a {N} criado')
            templates[N] = LDPC(dv, dc, N)

        num_N = len(possiveis_N)
        num_P = len(possiveis_P)
        monte_carlo_runs = 10000

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

# main()
# LDPC(3,7,98)
# print(Lista_N(dc=2))
# dc = 7 e dv=3

# [all_vnodes, all_cnodes]= LDPC(3,7,14)
# for i, vnode in enumerate(all_vnodes):
#     c_indices = [ all_cnodes.index(cnode) for cnode in vnode.c_nodes ]
#     print(f"VNode {i}: conectado a CNodes {c_indices}")
# print(all_vnodes)
# print(all_cnodes)


