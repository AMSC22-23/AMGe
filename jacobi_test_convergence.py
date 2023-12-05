# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def residual(phi, b):
    r = np.zeros_like(phi)
    for i in range(1, phi.shape[0] - 1):
        for j in range(1, phi.shape[1] - 1):
            r[i, j] = (
                4 * phi[i, j]
                - phi[i, j - 1]
                - phi[i - 1, j]
                - phi[i + 1, j]
                - phi[i, j + 1]
                - b[i, j]
            )
    return r

# ugly globals to keep count of total jacobi cost in multigrid
jacobi_iter = 0
jacobi_cost = 0

def jacobi(phi, b, iter, toll):
    global jacobi_iter, jacobi_cost
    ok = False
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
            jacobi_iter += k
            # cost of jacobi in phi.size**2 per iteration in theroy 
            # due to the matrix-vector product
            # since here it is matrix-free it has linear cost
            jacobi_cost += k * phi.size
            ok = True
            break
    if not ok:
        jacobi_iter += iter
        jacobi_cost += iter * phi.size
    return phi


def restrict(x):
    N = x.shape[0]
    n = (N - 1) // 2 + 1
    y = np.zeros((n, n))
    for i in range(2, N - 1, 2):
        for j in range(2, N - 1, 2):
            ii, jj = i // 2, j // 2
            y[ii, jj] = (
                0.25 * x[i, j]
                + 0.125 * (x[i + 1, j] + x[i, j + 1] + x[i - 1, j] + x[i, j - 1])
                + 0.0625
                * (
                    x[i + 1, j + 1]
                    + x[i - 1, j + 1]
                    + x[i - 1, j - 1]
                    + x[i + 1, j - 1]
                )
            )
    return y


def prolong(x):
    n = x.shape[0]
    N = (n - 1) * 2 + 1
    y = np.zeros((N, N))
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            if (i % 2 == 0) and (j % 2 == 0):
                y[i, j] = x[i // 2, j // 2]
            else:
                ii, jj = i // 2, j // 2
                y[i, j] = +0.25 * (
                    x[ii + 1, jj] + x[ii, jj + 1] + x[ii - 1, jj] + x[ii, jj - 1]
                ) + 0.125 * (
                    x[ii + 1, jj + 1]
                    + x[ii - 1, jj + 1]
                    + x[ii - 1, jj - 1]
                    + x[ii + 1, jj - 1]
                )
    return y


def mg_iter(phi, b, nu, toll, lvls):
    if lvls:
        phi = jacobi(phi, b, nu, toll)
        r = residual(phi, b)    
        rc = restrict(r)
        ec = np.zeros_like(rc)
        ec = mg_iter(ec, rc, nu, toll, lvls - 1)
        e = prolong(ec)
        phi = phi - e
        phi = jacobi(phi, b, nu, toll)
    else:
        phi = jacobi(phi, b, 100000, toll)
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

    phi0 = np.copy(phi)

    print("Building RHS")
    # RHS
    b = f(X, Y) * h * h

    toll = 1e-10
    for nlvls in range(int(np.log2(N)) - 1):
        print(f"Solving linear system with MG at {nlvls+1} levels")
        global jacobi_iter, jacobi_cost
        jacobi_iter = 0
        jacobi_cost = 0
        phi = np.copy(phi0)

        for _ in range(100000):
            phi = mg_iter(phi, b, 3, toll, nlvls)
            res = np.linalg.norm(residual(phi, b))
            if res < toll:
                break

        print(f"total iter {nlvls+1}lvl mg:", jacobi_iter)
        print(f"total cost {nlvls+1}lvl mg:", jacobi_cost)
        print("===================================")

    return X, Y, phi


def convergenge(Ns, u, f):
    e = []
    for N in Ns:
        X, Y, w = solve(N, u, f)
        U = u(X, Y)
        e.append(np.abs(w - U).max())
        print("error:", e[-1])
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
        return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0

    # force
    def f(x, y):
        return -(x * x * x + y * y * y)

    convergenge([32], u, f)
