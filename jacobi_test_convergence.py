# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def residual(phi, b):
    r = np.zeros_like(phi)
    for i in range(1, phi.shape[0] - 1):
        for j in range(1, phi.shape[1] - 1):
            r[i, j] =  (
                    4 * phi[i, j]
                    - phi[i, j - 1]
                    - phi[i - 1, j]
                    - phi[i + 1, j]
                    - phi[i, j + 1]
                    - b[i, j]
            )
    return r

def jacobi(phi, b, iter, toll):
    for k in range(iter):
        for i in range(1, phi.shape[0] - 1):
            for j in range(1, phi.shape[1] - 1):
                phi[i, j] = 0.25 * (
                    phi[i, j - 1]
                    + phi[i - 1, j]
                    + phi[i + 1, j]
                    + phi[i, j + 1]
                    + b[i, j]
                )
        r = residual(phi, b)
        err = np.linalg.norm(r)
        if err < toll:
            print(f"Solved in {k} iter")
            print("res:", err)
            break
    return phi


def solve(N, u, f):
    h = 1 / N
    x = np.linspace(0, 1.0, N + 1)
    y = np.linspace(0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y)

    # w is the full solution
    phi = np.zeros((N + 1, N + 1))

    # fix the boundary conditions
    phi[:, 0] = u(x, y[0])    # left Boundary
    phi[:, -1] = u(x, y[-1])  # Right Boundary
    phi[0, :] = u(x[0], y)    # Lower Boundary
    phi[-1, :] = u(x[-1], y)  # Upper Boundary

    print("Building RHS")
    # RHS
    b = f(X, Y) * h * h

    print("Solving linear system")
    jacobi(phi, b, 10000, 1e-14)

    return X, Y, phi


def convergenge(Ns, u, f):
    e = []
    for N in Ns:
        X, Y, w = solve(N, u, f)
        U = u(X, Y)
        e.append(np.abs(w - U).max())
        print("err:", e[-1])
        print("=================")
    plt.plot(Ns, e, "o-")
    # check that we converge with order 2 as by theory
    plt.plot(Ns, 5 / np.array(Ns) ** 2, "k--")
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("convergence_rate.png")


if __name__ == "__main__":
    # exact solution
    def u(x, y):
        return np.cos(2.0 * np.pi * x) * np.cos(2.0 * np.pi * y)

    # force
    def f(x, y):
        return 8.0 * np.pi * np.pi * np.cos(2.0 * np.pi * x) * np.cos(2.0 * np.pi * y)

    convergenge([4, 8, 16, 32], u, f)
