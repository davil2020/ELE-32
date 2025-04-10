import random

class VNode:
    def __init__(self, val: int = 0):
        if val not in (0, 1):
            raise ValueError("val deve ser 0 ou 1")
        self.val = val
        self.c_nodes = []
        self.dv = 0

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

def BSC_bit_flip(N, p):
    recieved = []
    for _ in range(N):
        recieved.append(1 if random.random() < p else 0) 
    return recieved


print(BSC_bit_flip(5,0.2))


# [all_vnodes, all_cnodes]= LDPC(3,7,7)
# print(all_vnodes)
# print(all_cnodes)


