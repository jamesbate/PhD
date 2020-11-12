import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def simplify_contours(qcset, thresh=plt.rcParams['path.simplify_threshold']):
    for c in qcset.collections:
        paths = c.get_paths()
        for i, p in enumerate(paths):
            p.simplify_threshold = thresh / plt.gcf().dpi
            paths[i] = p.cleaned(simplify=True)

# quadrupole tensor of positive charges in xy-direction and negative in -xy
Q = np.array([[0, 6., 0], [6., 0, 0], [0, 0, 0]])

xmax = 3.2
vmax = 0.15
ngrid = 400
levels = np.linspace(-vmax, vmax, 22)
XYZ = np.mgrid[-xmax:xmax:ngrid*1j, -xmax:xmax:ngrid*1j, 0:0:1j]
X, Y = XYZ[0,:,:,0], XYZ[1,:,:,0]

V = np.zeros_like(X)
for i in range(ngrid):
    for j in range(ngrid):
        r = XYZ[:,i,j,0]
        d = np.linalg.norm(r)
        V[i,j] = 0.5 * r.dot(Q).dot(r) / (4. * np.pi * d**5)
V = np.clip(V, -2.*vmax, 2.*vmax)

plt.figure(figsize=(6, 6)).add_axes([0, 0, 1, 1])
contf = plt.contourf(X, Y, V, levels=levels, cmap='RdBu_r', extend='both',
       norm=colors.SymLogNorm(linthresh=0.07, vmin=-vmax, vmax=vmax))
simplify_contours(contf, 0.9*plt.rcParams['lines.linewidth'])
cont = plt.contour(X, Y, V, levels=contf.levels, colors='k', linestyles='solid')
plt.xticks([]), plt.yticks([])
plt.gca().set_aspect(aspect='equal')
plt.gca().axis('off')
for i in -1,1:
    for j in -1,1:
        plt.text(.2*i, .2*j, {-1:u'\u2212', 1:'+'}[i*j],
                 size=18, ha='center', va='center')
plt.savefig('Quadrupole-potential-contour-xy.svg')