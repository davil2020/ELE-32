import numpy as np
import matplotlib as plt
import random

class VNode:
    def __init__(self):
        self.canal = 0
        self.c_nodes = []
        self.mensagens_entrada = []
        self.mensagens_saida = []
        self.dv = 0
        self.count_odd_c_node = 0

    def add_c_node(self, cnode):
        if not isinstance(cnode, CNode):
            raise TypeError("Esperado um objeto do tipo CNode")
        self.c_nodes.append(cnode)
        self.dv = len(self.c_nodes)


class CNode:
    def __init__(self):
        self.mensagens_saida = []
        self.mensagens_saida = []
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


num_bits = 1000
r = np.zeros(num_bits, dtype=int)

symbols = r + 1

SNR_db = 5
SNR_linear = 10**(SNR_db/10)

noise_variance = 1/ (2 * SNR_linear)
noise = np.sqrt(noise_variance) * np.random.randn(num_bits)

received = symbols + noise

#print(received[:20])

LLR = 2 * received / noise_variance

print(LLR[:100])
[v_nodes, c_nodes] = LDPC(3, 7, num_bits)

for i in range(1, 1000):
    v_nodes[i].canal = LLR[i] 

for v_node in v_nodes:
    for c_node in v_node.c_nodes:
        
