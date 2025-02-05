import numpy as np
import sys

a = np.array([1.0, 0.0, 0.0])
b = np.array([0.0, 1.0, 0.0])
c = np.array([0.0, 0.0, 1.0])

basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])

d_phi = 0.001
depth = 10

# print('symfind - simple crystal structure symmetry finder')

# if len(sys.argv) < 2:
#     print('usage: symfind vasp_file')
#     exit()

# with open(sys.argv[1], 'r') as file:
#     file

def rotate_x(A, phi):
    R_x = np.array([[1.0, 0.0,          0.0],
                    [0.0, np.cos(phi), -np.sin(phi)],
                    [0.0, np.sin(phi),  np.cos(phi)]])
    return np.dot(R_x, A)

def rotate_y(A, phi):
    R_x = np.array([[ np.cos(phi), 0.0, np.sin(phi)],
                    [ 0.0,         1.0, 0.0],
                    [-np.sin(phi), 0.0, np.cos(phi)]])
    return np.dot(R_x, A)

def rotate_z(A, phi):
    R_x = np.array([[np.cos(phi), -np.sin(phi), 0.0],
                    [np.sin(phi),  np.cos(phi), 0.0],
                    [0.0,          0.0,         np.cos(phi)]])
    return np.dot(R_x, A)

def to_cart(direct):
    return np.transpose(np.dot(np.transpose(np.array([a, b, c])), np.transpose(direct)))

def to_direct(cart):
    return np.transpose(np.dot(np.linalg.inv(np.transpose(np.array([a, b, c]))), np.transpose(cart)))

def reduce(direct):
    for v in direct:
        v[0] = v[0] % 1
        v[1] = v[1] % 1
        v[2] = v[2] % 1
        if v[0] < 0.0:
            v[0] += 1.0
        if v[1] < 0.0:
            v[1] += 1.0
        if v[2] < 0.0:
            v[2] += 1.0
    return direct

for i in range(depth):
    for j in range(depth):
        for k in range(depth):
            print('s')

print(to_cart(basis))
print(to_direct(to_cart(basis)))