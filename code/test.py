from qutip import sigmaz, tensor, qeye

z1 = tensor(qeye(2), sigmaz(), qeye(2))
z2 = tensor(sigmaz(), qeye(2), qeye(2))
z3 = tensor(qeye(2), qeye(2), sigmaz())

print(z1+z2 + z3)