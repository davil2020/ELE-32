class VNode:
    def __init__(self, val: int = 0):
        if val not in (0, 1):
            raise ValueError("val deve ser 0 ou 1")
        self.val = val
        self.v_nodes = []
        self.dv = 0

    def add_v_node(self, vnode):
        if not isinstance(vnode, VNode):
            raise TypeError("Esperado um objeto do tipo VNode")
        self.v_nodes.append(vnode)
        self.dv = len(self.v_nodes)


class CNode:
    def __init__(self):
        self.v_nodes = []
        self.dc = 0

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
            vnode.add_v_node(cnode)
            idx += 1

    return all_vnodes, all_cnodes

LDPC(3,7,7)
