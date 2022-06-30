"""
Composite graph utility
"""

# Modules
import pyamg
import networkx as nx
import numpy as np
import numpy.linalg as npla
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# Functions
def composite_graph(partition, amg, delta = 1):
    """
    Composite graph construction from AMG hierarchy
    - partition (numpy.array): Boolean array with 1s denoting the home nodes
    - amg (pyamg.levels): AMG hierarchy providing linear operators and integrid transfer operators
    - delta (int): Size of the expansion at every level
    """

    # AMG Hierarchy
    num_levels = len(amg.levels)
    num_nodes = [level.A.shape[0] for level in amg.levels]

    # Composite grid
    D = [np.zeros(num_nodes[l], dtype = int) for l in range(num_levels)]
    D[0][partition > 0] = 1.0

    for l in range(num_levels):
        A_l = amg.levels[l].A
        partition_tilde = D[l].copy()

        for nu in range(delta):
            partition_tilde[(A_l * partition_tilde != 0.0) | (partition_tilde != 0.0)] = 1

        if l == num_levels - 1:
            partition_tilde[:] = 1

        D[l][(D[l] == 0) & (partition_tilde > 0)] = 2

        if sum(D[l] == 0) == 0:
            break

        if l < num_levels - 1:
            P_l = amg.levels[l].P

            for row in range(P_l.shape[0]):
                if P_l.indptr[row + 1] - P_l.indptr[row] == 1:
                    if D[l][row] > 0:
                        col = P_l.indices[P_l.indptr[row]]
                        D[l + 1][col] = 1

    num_comp_levels = l + 1
    num_home = [sum(D[l] == 1) for l in range(num_comp_levels)]
    num_overlap = [sum(D[l] == 2) for l in range(num_comp_levels)]
    num_complement = [sum(D[l] == 0) for l in range(num_comp_levels)]
    num_comp_overlap = list(num_overlap)
    num_comp_overlap[0] += num_home[0]

    # Coarse to fine interpolator
    P_c = [None for _ in range(num_comp_levels - 1)]
    R_c = [None for _ in range(num_comp_levels - 1)]

    for l in range(num_comp_levels - 1, 0, - 1):
        P_l = amg.levels[l - 1].P

        # Mark fine nodes
        fine_nodes = np.zeros(num_nodes[l - 1], dtype = int)

        if l - 1 == 0:
            fine_nodes[D[l - 1] == 1] = np.arange(1, num_home[l - 1] + 1, dtype = int)
            fine_nodes[D[l - 1] == 2] = num_home[l - 1] + np.arange(1, num_overlap[l - 1] + 1, dtype = int)
            fine_nodes[D[l - 1] == 0] = num_home[l - 1] + num_overlap[l - 1] + np.arange(1, num_complement[l - 1] + 1, dtype = int)

        else:
            fine_nodes[D[l - 1] == 2] = np.arange(1, num_overlap[l - 1] + 1, dtype = int)
            fine_nodes[D[l - 1] == 0] = num_overlap[l - 1] + np.arange(1, num_complement[l - 1] + 1)

        # Mark coarse nodes
        coarse_nodes = np.zeros(num_nodes[l], dtype = int)

        if l - 1 == 0:
            coarse_nodes[(D[l] == 2) | (D[l] == 0)] = num_home[l - 1] + num_overlap[l - 1] + np.arange(1, num_overlap[l] + num_complement[l] + 1, dtype = int)
        else:
            coarse_nodes[(D[l] == 2) | (D[l] == 0)] = num_overlap[l - 1] + np.arange(1, num_overlap[l] + num_complement[l] + 1, dtype = int)

        for row in range(num_nodes[l - 1]):
            if P_l.indptr[row + 1] - P_l.indptr[row] == 1:
                ptr = P_l.indptr[row]
                col = P_l.indices[ptr]

                if l - 1 == 0:
                    if fine_nodes[row] <= (num_home[l - 1] + num_overlap[l - 1]):
                        coarse_nodes[col] = fine_nodes[row]

                else:
                    if fine_nodes[row] <= num_overlap[l - 1]:
                        coarse_nodes[col] = fine_nodes[row]

        # Construct level interpolator
        P_c_l_row = []
        P_c_l_col = []
        P_c_l_val = []

        if l - 1 == 0:
            num_fine = num_home[l - 1] + num_overlap[l - 1]
        else:
            num_fine = num_overlap[l - 1]

        for row in range(num_nodes[l - 1]):
            if fine_nodes[row] == 0:
                continue

            if fine_nodes[row] <= num_fine:
                P_c_l_row.append(fine_nodes[row] - 1)
                P_c_l_col.append(fine_nodes[row] - 1)
                P_c_l_val.append(1.0)

            else:
                for ptr in range(P_l.indptr[row], P_l.indptr[row + 1]):
                    col = P_l.indices[ptr]
                    val = P_l.data[ptr]

                    if coarse_nodes[col] > 0:
                        P_c_l_row.append(fine_nodes[row] - 1)
                        P_c_l_col.append(coarse_nodes[col] - 1)
                        P_c_l_val.append(val)

        if l - 1 == 0:
            P_c[l - 1] = spsp.csr_matrix((P_c_l_val, (P_c_l_row, P_c_l_col)), shape = (num_nodes[l - 1], num_home[l - 1] + num_overlap[l - 1] + num_overlap[l] + num_complement[l]))
        else:
            P_c[l - 1] = spsp.csr_matrix((P_c_l_val, (P_c_l_row, P_c_l_col)), shape = (num_overlap[l - 1] + num_complement[l - 1], num_overlap[l - 1] + num_overlap[l] + num_complement[l]))

        # Construct mapping to original ordering
        R_c_l_row = []
        R_c_l_col = []
        R_c_l_val = []
        dof = 0

        for row in range(num_nodes[l - 1]):
            if fine_nodes[row] > 0:
                R_c_l_row.append(dof)
                R_c_l_col.append(fine_nodes[row] - 1)
                R_c_l_val.append(1.0)
                dof += 1

        if l - 1 == 0:
            R_c[l - 1] = spsp.csr_matrix((R_c_l_val, (R_c_l_row, R_c_l_col)), shape = (num_nodes[l - 1], num_nodes[l - 1]))
        else:
            R_c[l - 1] = spsp.csr_matrix((R_c_l_val, (R_c_l_row, R_c_l_col)), shape = (num_overlap[l - 1] + num_complement[l - 1], num_overlap[l - 1] + num_complement[l - 1]))

    # Construct composite to global interpolator
    for l in range(num_comp_levels - 2, 0, - 1):
        I_c_lm1    = P_c[l - 1][:num_comp_overlap[l - 1], :num_comp_overlap[l - 1]]
        O_c_lm1    = P_c[l - 1][:num_comp_overlap[l - 1], num_comp_overlap[l - 1]:]
        P_c_lm1_21 = P_c[l - 1][num_comp_overlap[l - 1]:, :num_comp_overlap[l - 1]]
        P_c_lm1_22 = P_c[l - 1][num_comp_overlap[l - 1]:, num_comp_overlap[l - 1]:]

        P_c_lm1_22 = P_c_lm1_22 * (R_c[l] * P_c[l])
        P_c[l - 1] = spsp.bmat([[I_c_lm1, spsp.csr_matrix((I_c_lm1.shape[0], P_c_lm1_22.shape[1]))], [P_c_lm1_21, P_c_lm1_22]], format = "csr")

    Q_c = R_c[0] * P_c[0]

    # Mapping from composite to fine nodes
    nodes_to_fine = [np.arange(num_nodes[0], dtype = int)] + [np.empty(num_nodes[l], dtype = int) for l in range(1, num_comp_levels)]

    for l in range(num_comp_levels - 1):
        P_l = amg.levels[l].P

        for row in range(num_nodes[l]):
            if P_l.indptr[row + 1] - P_l.indptr[row] == 1:
                col = P_l.indices[P_l.indptr[row]]
                nodes_to_fine[l + 1][col] = nodes_to_fine[l][row]

    nodes_to_dofs = [np.zeros(num_nodes[l], dtype = int) for l in range(num_comp_levels)]
    nodes_to_dofs[0][D[0] == 1] = np.arange(1, num_home[0] + 1, dtype = int)
    nodes_to_dofs[0][D[0] == 2] = num_home[0] + np.arange(1, num_overlap[0] + 1, dtype = int)

    offset = num_home[0] + num_overlap[0]

    for l in range(num_comp_levels - 1):
        nodes_to_dofs[0][nodes_to_fine[l + 1][D[l + 1] == 2]] = offset + np.arange(1, num_overlap[l + 1] + 1, dtype = int)
        offset += num_overlap[l + 1]

    return Q_c, nodes_to_dofs[0]

# Testing
if __name__ == "__main__":
    # Create stiffness matrix
    m = 9
    n = m * m
    stencil = [[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]]
    A = pyamg.gallery.stencil_grid(stencil, (m, m), dtype = float, format = "csr")

    # Create AMG hierarchy
    amg = pyamg.ruge_stuben_solver(A)

    # Crate tagging
    partition = np.zeros((m, m) , dtype = int)
    partition[3:6, 3:6] = 1
    partition = partition.reshape(-1)

    # Compute composite assembly operator
    delta = 1
    Q_c, nodes_to_dofs = composite_graph(partition, amg, delta)

    # Compute composite operator
    A_c = Q_c.T * A * Q_c

    # Discretization grid
    h = 1.0 / (m + 1)
    x = np.kron(np.ones(m), np.linspace(h, 1.0 - h, m))
    y = np.kron(np.linspace(h, 1.0 - h, m), np.ones(m))

    R_hat = spsp.csr_matrix((np.ones(m), (np.arange(1, m + 1), np.arange(m))), shape = (m + 2, m))
    R = spsp.kron(R_hat, R_hat)

    # Draw graphs
    pos_fine = { i: (x[i], y[i]) for i in range(n) }
    pos_comp = { j - 1: (x[i], y[i]) for i, j in enumerate(nodes_to_dofs) if j > 0 }

    fig, ax = plt.subplots(figsize = (8, 8))
    nx.draw_networkx(nx.from_scipy_sparse_matrix(A), pos = pos_fine)
    ax.set_aspect("equal")
    ax.set_title("Fine operator graph")
    plt.savefig("fine_operator_graph.pdf", bbox_inches = "tight")

    fig, ax = plt.subplots(figsize = (8, 8))
    nx.draw_networkx(nx.from_scipy_sparse_matrix(A_c), pos = pos_comp)
    ax.set_aspect("equal")
    ax.set_title("Composite operator graph")
    plt.savefig("composite_operator_graph.pdf", bbox_inches = "tight")

    # Poisson problem
    f = 2.0 * (np.pi ** 2.0) * np.sin(np.pi * x) * np.sin(np.pi * y)
    u_star = np.sin(np.pi * x) * np.sin(np.pi * y)

    u = spspla.spsolve(A, (h ** 2.0) * f)
    v = Q_c * spspla.spsolve(A_c, Q_c.T * ((h ** 2.0) * f))

    X, Y = np.meshgrid(np.linspace(0.0, 1.0, m + 2), np.linspace(0.0, 1.0, m + 2))

    fig = plt.figure(figsize = (12, 8))
    ax = plt.axes(projection = "3d")
    ax.set_title("Exact solution")
    ax.plot_surface(X, Y, (R * u_star).reshape(m + 2, m + 2))
    plt.savefig("exact_solution.pdf", bbox_inches = "tight")

    fig = plt.figure(figsize = (12, 8))
    ax = plt.axes(projection = "3d")
    ax.set_title("Numerical solution with fine operator")
    ax.plot_surface(X, Y, (R * u).reshape(m + 2, m + 2))
    plt.savefig("numerical_solution_fine.pdf", bbox_inches = "tight")

    fig = plt.figure(figsize = (12, 8))
    ax = plt.axes(projection = "3d")
    ax.set_title("Numerical solution with composite operator")
    ax.plot_surface(X, Y, (R * v).reshape(m + 2, m + 2))
    plt.savefig("numerical_solution_composite.pdf", bbox_inches = "tight")

    plt.show()
