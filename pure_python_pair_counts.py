import numpy as np 

def python_npairs(x1, y1, z1, x2, y2, z2, rbins, Lbox):
    counts = np.zeros(len(rbins), dtype=int)
    rbins = rbins**2
    for i in range(len(x1)):
        for j in range(len(x2)):
            dx = np.minimum(np.fabs(x1[i] - x2[j]), 
                Lbox - np.fabs(x1[i] - x2[j]))
            dy = np.minimum(np.fabs(y1[i] - y2[j]), 
                Lbox - np.fabs(y1[i] - y2[j]))
            dz = np.minimum(np.fabs(z1[i] - z2[j]), 
                Lbox - np.fabs(z1[i] - z2[j]))
            dsq = dx*dx + dy*dy + dz*dz

            k = len(rbins)-1
            while dsq <= rbins[k]:
                counts[k] += 1
                k=k-1
                if k<0: break
    return counts

