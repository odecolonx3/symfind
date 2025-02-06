import numpy as np
import sys

a = np.array([1.0, 0.0, 0.0])
b = np.array([0.0, 1.0, 0.0])
c = np.array([0.0, 0.0, 1.0])

struct = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
basis   = np.array([a, b, c])

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

def J(A, basis):
    return np.transpose(np.dot(np.transpose(basis), np.transpose(A)))

def J_(A_, basis):
    return np.transpose(np.dot(np.linalg.inv(np.transpose(basis)), np.transpose(A_)))

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

def symcheck(struct, struct_):
    diff_tot = 0.0
    for i in range(len(struct)):
        diff = np.inf
        for j in range(len(struct_)):
            if (np.linalg.norm(struct[i] - struct_[j]) < diff):
                diff = np.linalg.norm(struct[i] - struct_[j])
            # print('symfind: symcheck(): ', i, ' --- ', j, ', diff = ', np.linalg.norm(struct[i] - struct_[j]))
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

            if ((i, j, k) != (1, 1, 1)) and ((i, j, k) != (1, 0, 0)):
                continue

            axis  = a * i + b * j + c * k
            beta  = np.arccos(np.dot(axis, np.array([0.0, 0.0, 1.0])) / (np.linalg.norm(axis)))
            alpha = np.arccos(axis[0] / (np.linalg.norm(axis) * np.cos(np.pi / 2 - beta)))

            # struct in basis -> struct in cart

            struct_cart = J(struct, basis)

            # basis -> basis_Rzx

            basis_Rz  = rotate_z(basis, -alpha)
            basis_Rzx = rotate_x(basis, -beta)

            # struct in cart -> struct in basis_Rzx

            struct_Rzx = J_(struct_cart, basis_Rzx)

            print('basis       = ', basis)
            print('struct      = ', struct)
            print('struct_cart = ', struct_cart)
            print('basis_Rzx   = ', basis_Rzx)
            print('struct_Rzx  = ', struct_Rzx)

            # struct_Rzx = to_cart(reduce(struct))
            # struct_Rzx = rotate_z(struct_Rzx, alpha)
            # struct_Rzx = rotate_x(struct_Rzx, beta)
            # struct_Rzx = to_direct(struct_Rzx)
            # struct_Rzx = reduce(struct_Rzx)
            # struct_Rzx = to_cart(struct_Rzx)

            # for ord in axis_order:
            #     gamma      = 2 * np.pi / ord
                # struct_Rzxz = rotate_z(struct_Rzx, gamma)

                # struct_Rzxz = to_direct(struct_Rzxz)
                # struct_Rzxz = reduce(struct_Rzxz)
                # struct_Rzxz = to_cart(struct_Rzxz)

                # if ((i, j, k, ord) == (1, 1, 1, 1)):
                #     print('struct      = ', to_cart(reduce(struct)))
                #     print('struct_Rzx  = ', struct_Rzx)
                #     print('struct_Rzxz = ', struct_Rzxz)
                #     print('struct_Rzxz[1][1] = ', struct_Rzxz[1][1] < 1.0)
                #     print('basis_Rzxz = ', rotate_z(rotate_x(rotate_z(np.array([a, b, c]), alpha), beta), gamma))

                # if (symcheck(struct_Rzx, struct_Rzxz)):
                #     print('symfind: found axial symmetry along [%d%d%d] axis of %d order' % (i, j, k, ord))



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

                        # struct_ = to_cart(reduce(struct))
                        # struct_ = rotate_x(struct_, phi_x)
                        # # struct_ = rotate_y(struct_, phi_y)
                        # # struct_ = rotate_z(struct_, phi_z)
                        # struct_ = to_direct(struct_)
                        # struct_ = reduce(struct_)
                        # struct_ = to_cart(struct_)
                        
#                         norm = 0.0

#                         for l in range(len(struct)):
#                             norm += np.sqrt((struct_[l][0] - to_cart(reduce(struct))[l][0])**2 + (struct_[l][1] - to_cart(reduce(struct))[l][1])**2 + (struct_[l][2] - to_cart(reduce(struct))[l][2])**2)

#                         if (norm < epsilon):
#                             print('symfind: point irreducible symmmetry rotation along OX found, phi = 2 pi / ', o_x, ' * ', i, ', order = ', o_x, ', with precision ~', epsilon)