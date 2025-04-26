import random

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

    connections = []

    for i in range(N):
        connections += [i] * dv

    import random
    random.shuffle(connections)

    idx = 0
    for cnode in all_cnodes:
        for _ in range(dc):
            vnode_index = connections[idx]
            vnode = all_vnodes[vnode_index]
            cnode.add_v_node(vnode)
            vnode.add_c_node(cnode)
            idx += 1

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

# print(Lista_N())

def LDPC_decode():
    dc = 7
    dv = 3
    possiveis_N = Lista_N(dc)

    N = possiveis_N[0]
    [all_vnodes, all_cnodes]= LDPC(dv, dc, N)
    BSC_bit_flip(all_vnodes, 0.1)

    # print([vnode.val for vnode in all_vnodes])

    return 0

LDPC_decode()

# dc = 7 e dv=3

# [all_vnodes, all_cnodes]= LDPC(3,7,7)
# print(all_vnodes)
# print(all_cnodes)


