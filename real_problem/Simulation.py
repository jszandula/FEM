import numpy as np

class GaussElimination:
    def gaussElimination(self, H, P):
        n = H.shape[0]
        
        for column in range(0, n-1): 
            for row in range(column+1, n): 
                l = (H[row, column] / H[column, column])
                H[row] = H[row] - H[column] * l
                P[row] = P[row] - P[column] * l

        x = np.zeros(n)
        x[n-1] = P.item(n-1) / H.item(n-1,n-1)
        for k in range(n-2, -1, -1):
            x[k] = (P[k] - sum(H[k] * np.matrix(x).T)).item(0) / H.item(k,k) 
        return x
      