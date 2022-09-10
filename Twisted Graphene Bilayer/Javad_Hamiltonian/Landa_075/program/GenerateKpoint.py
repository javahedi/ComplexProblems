import numpy as np

def trisample(A, B, C):
    """
    Given three vertices A, B, C, 
    sample point uniformly in the triangle
    """
    r1 = np.random.random()
    r2 = np.random.random()

    s1 = np.sqrt(r1)

    x = A[0] * (1.0 - s1) + B[0] * (1.0 - r2) * s1 + C[0] * r2 * s1
    y = A[1] * (1.0 - s1) + B[1] * (1.0 - r2) * s1 + C[1] * r2 * s1

    return (x, y)
