import numpy as np
import sys
import matplotlib.pyplot as plt

a = np.array([1.0, 0.0, 0.0])
b = np.array([0.0, 1.0, 0.0])
c = np.array([0.0, 0.0, 1.0])

struct = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
basis   = np.array([a, b, c])

epsilon = 0.001

axis_order = [1, 2, 3, 4, 6]

depth = 4
scale = 2

print('symfind - simple crystal structure symmetry finder')

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
    R_y = np.array([[ np.cos(phi), 0.0, np.sin(phi)],
                    [ 0.0,         1.0, 0.0],
                    [-np.sin(phi), 0.0, np.cos(phi)]])
    return np.transpose(np.dot(R_y, np.transpose(A)))

def rotate_z(A, phi):
    R_z = np.array([[np.cos(phi), -np.sin(phi), 0.0],
                    [np.sin(phi),  np.cos(phi), 0.0],
                    [0.0,          0.0,         1.0]])
    return np.transpose(np.dot(R_z, np.transpose(A)))

def J(A, basis):
    return np.transpose(np.dot(np.transpose(basis), np.transpose(A)))

def J_(A_, basis):
    return np.transpose(np.dot(np.linalg.inv(np.transpose(basis)), np.transpose(A_)))

def reduce(struct):
    for v in struct:
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

    return struct

def extend(struct):
    tmp = np.zeros((len(struct) * 8, 3))
    for j in range(0, 8):
        for i in range(len(struct) * j, len(struct) * (j + 1)):
            if (j == 0):
                tmp[i] = struct[i % len(struct)]
            if (j == 1):
                tmp[i] = struct[i % len(struct)] - [1.0, 0.0, 0.0]
            if (j == 2):
                tmp[i] = struct[i % len(struct)] - [0.0, 1.0, 0.0]
            if (j == 3):
                tmp[i] = struct[i % len(struct)] - [0.0, 0.0, 1.0]
            if (j == 4):
                tmp[i] = struct[i % len(struct)] - [1.0, 1.0, 0.0]
            if (j == 5):
                tmp[i] = struct[i % len(struct)] - [0.0, 1.0, 1.0]
            if (j == 6):
                tmp[i] = struct[i % len(struct)] - [1.0, 0.0, 1.0]
            if (j == 7):
                tmp[i] = struct[i % len(struct)] - [1.0, 1.0, 1.0]
    return tmp

def symcheck(struct, struct_, basis, basis_):
    diff_tot = 0.0
    
    struct_basis  = J_(struct,  basis)
    struct_basis_ = J_(struct_, basis_)

    for i in range(len(struct)):
        if (np.abs(struct_basis[i][0]) >= 1.0) or \
            (np.abs(struct_basis[i][1]) >= 1.0) or \
                (np.abs(struct_basis[i][2]) >= 1.0) or \
                    (np.abs(np.abs(struct_basis[i][0]) - 1.0) < epsilon) or \
                        (np.abs(np.abs(struct_basis[i][1]) - 1.0) < epsilon) or \
                            (np.abs(np.abs(struct_basis[i][2]) - 1.0) < epsilon):
            continue

        diff = np.inf

        for j in range(len(struct_)):

            if (np.abs(struct_basis_[j][0]) >= 1.0) or \
                (np.abs(struct_basis_[j][1]) >= 1.0) or \
                    (np.abs(struct_basis_[j][2]) >= 1.0) or \
                        (np.abs(np.abs(struct_basis_[j][0]) - 1.0) < epsilon) or \
                            (np.abs(np.abs(struct_basis_[j][1]) - 1.0) < epsilon) or \
                                (np.abs(np.abs(struct_basis_[j][2]) - 1.0) < epsilon):
                continue

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
        # print('symfind: symcheck(): diff_tot = %.8f < epsilon = %.8f, symmetry check passed!' % (diff_tot, epsilon))
        return True
    else:
        # print('symfind: symcheck(): diff_tot = %.8f > epsilon = %.8f, symmetry check failed...' % (diff_tot, epsilon))
        return False

# Reducing and extending the structure

struct = reduce(struct)
struct = extend(struct)

for i in range(depth):
    for j in range(depth):
        for k in range(depth):
            if ((i == 0) and (j == 0) and (k == 0)):
                continue
            if (np.gcd.reduce([i, j, k]) != 1):
                continue

            print('\nsymfind: checking [%d%d%d] axis... \n' % (i, j, k))

            # if ((i, j, k) != (1, 0, 2)):
            #     continue

            # rotational axis, defined as Weiss indexes

            axis = a * i + b * j + c * k
        
            # axis in basis -> axis in cart

            axis_cart = J(axis, basis)

            # Euler precession and nutation angles

            beta  = np.arccos(np.dot(axis_cart, np.array([0.0, 0.0, 1.0])) / (np.linalg.norm(axis_cart)))
            alpha = np.arccos(axis_cart[0] / (np.linalg.norm(axis_cart) * np.cos(np.pi / 2 - beta)))

            # print('alpha = %.2f pi' % (alpha / np.pi))
            # print('beta  = %.2f pi' % (beta  / np.pi))

            # axis in cart -> axis (rotated Rzx) in cart

            axis_cart_Rz  = rotate_z(axis_cart, alpha)
            axis_cart_Rzx = rotate_x(axis_cart_Rz, beta)

            # struct in basis -> struct in cart

            struct_cart = J(struct, basis)

            # struct in cart -> struct (rotated Rzx) in cart 

            struct_cart_Rz  = rotate_z(struct_cart, alpha)
            struct_cart_Rzx = rotate_x(struct_cart_Rz, beta)
           
            # basis -> basis (rotated Rzx)

            basis_Rz  = rotate_z(basis, alpha)
            basis_Rzx = rotate_x(basis_Rz, beta)

            # struct (rotated Rzx) in cart -> struct (rotated Rzx) in basis (rotated Rzx)

            struct_Rzx = J_(struct_cart_Rzx, basis_Rzx)

            # print('basis            = ', basis)
            # print('basis_Rz         = ', basis_Rz)
            # print('basis_Rzx        = ', basis_Rzx)
            # print('axis             = ', axis_cart)
            # print('axis_cart        = ', axis_cart)
            # print('axis_cart_Rz     = ', axis_cart_Rz)
            # print('axis_cart_Rzx    = ', axis_cart_Rzx)
            # print('struct           = ', struct)
            # print('struct_cart      = ', struct_cart)
            # print('struct_cart_Rz   = ', struct_cart_Rz)
            # print('struct_cart_Rzx  = ', struct_cart_Rzx)
            # print('struct_Rzx       = ', struct_Rzx)

            # Checking symmetries of different orders

            for ord in axis_order:

                # if ord != 3:
                #     continue

                # Euler rotation angle

                gamma = 2 * np.pi / ord

                # struct (rotated Rzx, reducted and extended) in cart -> struct (rotated Rzxz, reducted and extended) in cart 

                struct_cart_Rzxz = rotate_z(struct_cart_Rzx, gamma)

                # basis (rotated Rzx) -> basis (rotated Rzxz)

                basis_Rzxz = rotate_z(basis_Rzx, gamma)

                # struct (rotated Rzxz, reducted and extended) in cart  -> struct (rotated Rzxz, reducted and extended) in basis (rotated Rzxz)

                struct_Rzxz = J_(struct_cart_Rzxz, basis_Rzxz)

                # print('basis_Rzxz        = ', basis_Rzxz)
                # print('struct_cart_Rzxz  = ', struct_cart_Rzxz)
                # print('struct_Rzxz       = ', struct_Rzxz)

                # fig = plt.figure()
                # ax = fig.add_subplot(projection='3d')
                # ax.scatter(struct_cart_Rzx[:, 0], struct_cart_Rzx[:, 1], struct_cart_Rzx[:, 2], marker='o')
                # ax.scatter(struct_cart_Rzxz[:, 0], struct_cart_Rzxz[:, 1], struct_cart_Rzxz[:, 2], marker='^')

                # ax.set_xlabel('X Label')
                # ax.set_ylabel('Y Label')
                # ax.set_zlabel('Z Label')

                # plt.show()

                if (symcheck(struct_cart_Rzx, struct_cart_Rzxz, basis_Rzx, basis_Rzxz)):
                    print('symfind: found point symmetry along [%d%d%d] axis of %d order!' % (i, j, k, ord))
                else:
                    print('symfind: %d order symcheck failed...' % (ord))