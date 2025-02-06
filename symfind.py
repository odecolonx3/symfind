import numpy as np
import sys

a = np.array([1.0, 0.0, 0.0])
b = np.array([0.0, 1.0, 0.0])
c = np.array([0.0, 0.0, 1.0])

basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])

epsilon = 0.01

axis_order = [1, 2, 3, 4, 6]

depth = 3
scale = 2

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
    return np.transpose(np.dot(R_x, np.transpose(A)))

def rotate_y(A, phi):
    R_x = np.array([[ np.cos(phi), 0.0, np.sin(phi)],
                    [ 0.0,         1.0, 0.0],
                    [-np.sin(phi), 0.0, np.cos(phi)]])
    return np.transpose(np.dot(R_x, np.transpose(A)))

def rotate_z(A, phi):
    R_x = np.array([[np.cos(phi), -np.sin(phi), 0.0],
                    [np.sin(phi),  np.cos(phi), 0.0],
                    [0.0,          0.0,         1.0]])
    return np.transpose(np.dot(R_x, np.transpose(A)))

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

        if v[0] > 1.0:
            v[0] -= 1.0
        if v[1] > 1.0:
            v[1] -= 1.0
        if v[2] > 1.0:
            v[2] -= 1.0

        if v[0] == 1.0:
            v[0] = 0.0
        if v[1] == 1.0:
            v[1] = 0.0
        if v[2] == 1.0:
            v[2] = 0.0
        
        if np.abs(v[0] - 1.0) < epsilon:
            v[0] = 0.0
        if np.abs(v[1] - 1.0) < epsilon:
            v[1] = 0.0
        if np.abs(v[2] - 1.0) < epsilon:
            v[2] = 0.0

    return direct

def expand(direct):
   print() 

def symcheck(basis, basis_):
    diff_tot = 0.0
    for i in range(len(basis)):
        diff = np.inf
        for j in range(len(basis_)):
            if (np.linalg.norm(basis[i] - basis_[j]) < diff):
                diff = np.linalg.norm(basis[i] - basis_[j])
            # print('symfind: symcheck(): ', i, ' --- ', j, ', diff = ', np.linalg.norm(basis[i] - basis_[j]))
        if (np.isinf(diff)):
            print('symfind: symcheck(): error! infinite norm')
            exit(1)
        # print('symfind: symcheck(): ', i, ', diff = ', diff)
        diff_tot += diff ** 2
    diff_tot = np.sqrt(diff_tot)
    if (diff_tot < epsilon):
        # print('symfind: symcheck(): symmetry check passed!')
        return True
    else:
        # print('symfind: symcheck(): symmetry check failed...')
        return False

for i in range(depth):
    for j in range(depth):
        for k in range(depth):
            if  ((i == 0) and (j == 0) and (k == 0)):
                continue
            if (np.gcd.reduce([i, j, k]) != 1):
                continue

            print('\nsymfind: checking [%d%d%d] axis... \n' % (i, j, k))

            axis  = a * i + b * j + c * k
            beta  = np.arccos(np.dot(axis, np.array([0.0, 0.0, 1.0])) / (np.linalg.norm(axis)))
            alpha = np.arccos(axis[0] / (np.linalg.norm(axis) * np.cos(np.pi / 2 - beta)))

            basis_Rzx = to_cart(reduce(basis))
            basis_Rzx = rotate_z(basis_Rzx, alpha)
            basis_Rzx = rotate_x(basis_Rzx, beta)
            basis_Rzx = to_direct(basis_Rzx)
            basis_Rzx = reduce(basis_Rzx)
            basis_Rzx = to_cart(basis_Rzx)

            for ord in axis_order:
                gamma      = 2 * np.pi / ord
                basis_Rzxz = rotate_z(basis_Rzx, gamma)

                basis_Rzxz = to_direct(basis_Rzxz)
                basis_Rzxz = reduce(basis_Rzxz)
                basis_Rzxz = to_cart(basis_Rzxz)

                if ((i, j, k, ord) == (1, 1, 1, 1)):
                    print('basis      = ', to_cart(reduce(basis)))
                    print('basis_Rzx  = ', basis_Rzx)
                    print('basis_Rzxz = ', basis_Rzxz)
                    print('basis_Rzxz[1][1] = ', basis_Rzxz[1][1] < 1.0)
                    print('abc_Rzxz = ', rotate_z(rotate_x(rotate_z(np.array([a, b, c]), alpha), beta), gamma))

                if (symcheck(basis_Rzx, basis_Rzxz)):
                    print('symfind: found axial symmetry along [%d%d%d] axis of %d order' % (i, j, k, ord))



# for o_x in axis_order:
#     # for o_y in axis_order:
#     #     for o_z in axis_order:
#             for i in range(o_x):
#                 # for j in range(o_y):
#                 #     for k in range(o_z):
#                         phi_x = 2 * np.pi / o_x * i
#                         # phi_y = 2 * np.pi / o_y * j
#                         # phi_z = 2 * np.pi / o_z * k

#                         if (i == 0) or ((o_x % i == 0) and (i != 1)):
#                             continue

                        # basis_ = to_cart(reduce(basis))
                        # basis_ = rotate_x(basis_, phi_x)
                        # # basis_ = rotate_y(basis_, phi_y)
                        # # basis_ = rotate_z(basis_, phi_z)
                        # basis_ = to_direct(basis_)
                        # basis_ = reduce(basis_)
                        # basis_ = to_cart(basis_)
                        
#                         norm = 0.0

#                         for l in range(len(basis)):
#                             norm += np.sqrt((basis_[l][0] - to_cart(reduce(basis))[l][0])**2 + (basis_[l][1] - to_cart(reduce(basis))[l][1])**2 + (basis_[l][2] - to_cart(reduce(basis))[l][2])**2)

#                         if (norm < epsilon):
#                             print('symfind: point irreducible symmmetry rotation along OX found, phi = 2 pi / ', o_x, ' * ', i, ', order = ', o_x, ', with precision ~', epsilon)