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
