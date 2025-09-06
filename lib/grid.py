import numpy as np
from scipy.sparse import coo_matrix
from matplotlib.tri import Triangulation


class Grid:
    """
    Encapsulates mesh generation, triangulation, FEM matrix assembly, and source vector setup.
    """

    def __init__(
        self,
        m: int,
        n: int,
        left: float = 0.0,
        right: float = 1.0,
        bottom: float = 0.0,
        top: float = 1.0,
    ):
        self.m = m
        self.n = n
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

        ## FIXME: properties here need to be adjusted based on changes to the corresponding methods
        self.tri = None
        self.A = None
        self.B = None
        self.RHS = None

    @property
    def get_spacing(self):
        return self.dx, self.dy

    @property
    def get_meshgrid(self):
        self.x, self.y = self.set_FD_meshgrid()
        return self.x, self.y

    @property
    def dx(self):
        return (self.right - self.left) / self.m

    @property
    def dy(self):
        return (self.top - self.bottom) / self.n

    def set_FD_meshgrid(self):
        """
        Generates FD coordinate grid for SP-flooding.
        Used for the transport equations.

        Returns:
            (x, y) meshgrid arrays
        """
        x = np.linspace(self.left, self.right, self.m + 1)
        y = np.linspace(self.bottom, self.top, self.n + 1)
        self.x, self.y = np.meshgrid(x, y)
        return self.x, self.y

    def get_flat_index_matrix(self) -> np.ndarray:
        """
        Returns a matrix of shape (n+1, m+1) with flat indices at each grid point.
        """
        return np.arange((self.m + 1) * (self.n + 1)).reshape((self.n + 1, self.m + 1))


class FEMesh(Grid):
    def __init__(self, m: int, n: int):
        super().__init__(m, n)
        self.U = None
        self.L = None
        self.grid_size = None
        self.right_hand = None
        self.A = None
        self.B = None
        self.sparsed_A = None

    def set_triangulation(self):
        #  Setting up triangulations for the FEM grid
        #  U = cell array with each element = array of vertices of Upper Triangle of
        #  the rectangular cell
        #  L = cell array with each element = array of vertices of Lower Triangle of
        #  the rectangular cell
        #  At every point (i,j), U{i,j} & L{i,j} are cells with coordinates of vertices
        #  of the two triangles obtained by bisecting the rectangle starting at
        #  (i,j). The bisection line goes from NW to SE.
        self.U = np.empty((self.m, self.n), dtype=object)
        self.L = np.empty((self.m, self.n), dtype=object)

        for j in range(self.m):
            for k in range(self.n):
                x1 = self.left + j * self.dx
                y1 = self.bottom + k * self.dy
                x2 = self.left + (j + 1) * self.dx
                y2 = y1
                x3 = x1
                y3 = self.bottom + (k + 1) * self.dy
                x4 = x2
                y4 = y3

                # lower triangle vertices
                l = {"x": np.array([x1, x2, x3]), "y": np.array([y1, y2, y3])}

                # upper triangle vertices
                u = {"x": np.array([x4, x3, x2]), "y": np.array([y4, y3, y2])}

                self.U[j, k] = u
                self.L[j, k] = l

    def _polyarea(self, x, y):
        """
        Calculate the area of a polygon using the Shoelace formula.
        The vertices are defined by the x and y coordinates.

        Parameters:
        x (list or array): x-coordinates of the polygon vertices
        y (list or array): y-coordinates of the polygon vertices

        Returns:
        float: Area of the polygon
        """
        return 0.5 * abs(
            sum(x[i] * y[i + 1] - y[i] * x[i + 1] for i in range(-1, len(x) - 1))
        )

    def _beta_func(self, x, y, beta):
        """
        Evaluates the coefficient beta at each grid point
        $$ \beta = K \lambda $$
        x and y are the coordinates of the grid point
        The corresponding index locations in the matrix for beta
        are determined in mm and nn respectively.

        Analogous to beta_func.m in the MATLAB code
        """
        mm = round((x - self.left) / self.dx)
        nn = round((y - self.bottom) / self.dy)

        return beta[
            int(nn), int(mm)
        ]  # nn and mm are casted to integers for indexing purposes

    def _set_FE_meshgrid_helper(self, T, beta, V):
        """
        Evaluates beta at the vertices of the element triangle

        Input:
        % T is a structure array with fields x & y where
        %   T.x contains x coordinates of vertices of an element triangle
        %   T.y contains y coordinates of vertices of an element triangle
        % beta is the average of the value at the vertices of the
        %   coefficient $$\beta = K(x) \lambda(s,c,\Gamma)$$

        Analogous to the weak.m function in the MATLAB code
        """
        beta_1 = self._beta_func(T["x"][0], T["y"][0], beta)
        beta_2 = self._beta_func(T["x"][1], T["y"][1], beta)
        beta_3 = self._beta_func(T["x"][2], T["y"][2], beta)

        # computing average of the beta values at the vertices
        beta_avg = (beta_1 + beta_2 + beta_3) / 3

        s = self._polyarea(T["x"], T["y"])

        # Create and manipulate matrix M
        M = np.vstack((T["x"], T["y"], [1, 1, 1])).T
        M_inv = np.linalg.inv(M)
        M = M_inv[:2, :]  # Extract the first two rows of M_inv

        # Calculate vdiff and inte
        vdiff = np.dot(M, V)
        inte = np.dot(vdiff.T, beta_avg * np.dot(M, s))

        # Output result
        inte = np.append(inte, [0])

        return inte

    def set_FE_meshgrid(self, beta):
        """
        Generate FE coordinate grid for elliptic pressure calculations
        Analogous to the setGrid.m function in the MATLAB code
        """
        self.grid_size = np.empty((self.m + 1, self.n + 1), dtype=object)

        for j in range(self.m + 1):
            for l in range(self.n + 1):

                if j == 0 and l != 0 and l != self.n:
                    t1 = self._set_FE_meshgrid_helper(
                        self.L[j, l], beta, np.array([1, 0, 0])
                    )
                    t2 = np.array([0, 0, 0, 0])
                    t3 = np.array([0, 0, 0, 0])
                    t4 = np.array([0, 0, 0, 0])
                    t5 = self._set_FE_meshgrid_helper(
                        self.L[j, l - 1], beta, np.array([0, 0, 1])
                    )
                    t6 = self._set_FE_meshgrid_helper(
                        self.U[j, l - 1], beta, np.array([0, 1, 0])
                    )

                if j == self.m and l != 0 and l != self.n:
                    t1 = np.array([0, 0, 0, 0])
                    t2 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l], beta, np.array([0, 0, 1])
                    )
                    t3 = self._set_FE_meshgrid_helper(
                        self.L[j - 1, l], beta, np.array([0, 1, 0])
                    )
                    t4 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l - 1], beta, np.array([1, 0, 0])
                    )
                    t5 = np.array([0, 0, 0, 0])
                    t6 = np.array([0, 0, 0, 0])

                if j != 0 and j != self.m and l == 0:
                    t1 = self._set_FE_meshgrid_helper(
                        self.L[j, l], beta, np.array([1, 0, 0])
                    )
                    t2 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l], beta, np.array([0, 0, 1])
                    )
                    t3 = self._set_FE_meshgrid_helper(
                        self.L[j - 1, l], beta, np.array([0, 1, 0])
                    )
                    t4 = np.array([0, 0, 0, 0])
                    t5 = np.array([0, 0, 0, 0])
                    t6 = np.array([0, 0, 0, 0])

                if j != 0 and j != self.m and l == self.n:
                    t1 = np.array([0, 0, 0, 0])
                    t2 = np.array([0, 0, 0, 0])
                    t3 = np.array([0, 0, 0, 0])
                    t4 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l - 1], beta, np.array([1, 0, 0])
                    )
                    t5 = self._set_FE_meshgrid_helper(
                        self.L[j, l - 1], beta, np.array([0, 0, 1])
                    )
                    t6 = self._set_FE_meshgrid_helper(
                        self.U[j, l - 1], beta, np.array([0, 1, 0])
                    )

                if j == 0 and l == 0:
                    t1 = self._set_FE_meshgrid_helper(
                        self.L[j, l], beta, np.array([1, 0, 0])
                    )
                    t2 = np.array([0, 0, 0, 0])
                    t3 = np.array([0, 0, 0, 0])
                    t4 = np.array([0, 0, 0, 0])
                    t5 = np.array([0, 0, 0, 0])
                    t6 = np.array([0, 0, 0, 0])

                if j == 0 and l == self.n:
                    t1 = np.array([0, 0, 0, 0])
                    t2 = np.array([0, 0, 0, 0])
                    t3 = np.array([0, 0, 0, 0])
                    t4 = np.array([0, 0, 0, 0])
                    t5 = self._set_FE_meshgrid_helper(
                        self.L[j, l - 1], beta, np.array([0, 0, 1])
                    )
                    t6 = self._set_FE_meshgrid_helper(
                        self.U[j, l - 1], beta, np.array([0, 1, 0])
                    )

                if j == self.m and l == 0:
                    t1 = np.array([0, 0, 0, 0])
                    t2 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l], beta, np.array([0, 0, 1])
                    )
                    t3 = self._set_FE_meshgrid_helper(
                        self.L[j - 1, l], beta, np.array([0, 1, 0])
                    )
                    t4 = np.array([0, 0, 0, 0])
                    t5 = np.array([0, 0, 0, 0])
                    t6 = np.array([0, 0, 0, 0])

                if j == self.m and l == self.n:
                    t1 = np.array([0, 0, 0, 0])
                    t2 = np.array([0, 0, 0, 0])
                    t3 = np.array([0, 0, 0, 0])
                    t4 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l - 1], beta, np.array([1, 0, 0])
                    )
                    t5 = np.array([0, 0, 0, 0])
                    t6 = np.array([0, 0, 0, 0])

                if j != 0 and j != self.m and l != 0 and l != self.n:
                    t1 = self._set_FE_meshgrid_helper(
                        self.L[j, l], beta, np.array([1, 0, 0])
                    )
                    t2 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l], beta, np.array([0, 0, 1])
                    )
                    t3 = self._set_FE_meshgrid_helper(
                        self.L[j - 1, l], beta, np.array([0, 1, 0])
                    )
                    t4 = self._set_FE_meshgrid_helper(
                        self.U[j - 1, l - 1], beta, np.array([1, 0, 0])
                    )
                    t5 = self._set_FE_meshgrid_helper(
                        self.L[j, l - 1], beta, np.array([0, 0, 1])
                    )
                    t6 = self._set_FE_meshgrid_helper(
                        self.U[j, l - 1], beta, np.array([0, 1, 0])
                    )

                # formulating grid
                grid = {
                    "c": t1[0] + t2[2] + t3[1] + t4[0] + t5[2] + t6[1],
                    "w": t3[0] + t4[1],
                    "s": t4[2] + t5[0],
                    "n": t1[2] + t2[0],
                    "e": t1[1] + t6[0],
                    "nw": t2[1] + t3[2],
                    "se": t5[1] + t6[2],
                    "const": t1[3] + t2[3] + t3[3] + t4[3] + t5[3] + t6[3],
                }

                self.grid_size[j, l] = grid

    def set_right_hand(self, source_prod_matrix):
        self.right_hand = np.zeros(((self.m + 1) * (self.n + 1), 1))

        for j in range(self.m + 1):
            for l in range(self.n + 1):

                # finding corresponding index
                idx = j + l * (self.m + 1)

                if j == 0 and l != 0 and l != self.n:
                    t1 = self._FInt(
                        self.L[j, l], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t2 = 0
                    t3 = 0
                    t4 = 0
                    t5 = self._FInt(
                        self.L[j, l - 1], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t6 = self._FInt(
                        self.U[j, l - 1], source_prod_matrix, np.array([0, 1, 0])
                    )

                if j == self.m and l != 0 and l != self.n:
                    t1 = 0
                    t2 = self._FInt(
                        self.U[j - 1, l], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t3 = self._FInt(
                        self.L[j - 1, l], source_prod_matrix, np.array([0, 1, 0])
                    )
                    t4 = self._FInt(
                        self.U[j - 1, l - 1], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t5 = 0
                    t6 = 0

                if j != 0 and j != self.m and l == 0:
                    t1 = self._FInt(
                        self.L[j, l], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t2 = self._FInt(
                        self.U[j - 1, l], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t3 = self._FInt(
                        self.L[j - 1, l], source_prod_matrix, np.array([0, 1, 0])
                    )
                    t4 = 0
                    t5 = 0
                    t6 = 0

                if j != 0 and j != self.m and l == self.n:
                    t1 = 0
                    t2 = 0
                    t3 = 0
                    t4 = self._FInt(
                        self.U[j - 1, l - 1], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t5 = self._FInt(
                        self.L[j, l - 1], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t6 = self._FInt(
                        self.U[j, l - 1], source_prod_matrix, np.array([0, 1, 0])
                    )

                if j == 0 and l == 0:
                    t1 = self._FInt(
                        self.L[j, l], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t2 = 0
                    t3 = 0
                    t4 = 0
                    t5 = 0
                    t6 = 0

                if j == 0 and l == self.n:
                    t1 = 0
                    t2 = 0
                    t3 = 0
                    t4 = 0
                    t5 = self._FInt(
                        self.L[j, l - 1], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t6 = self._FInt(
                        self.U[j, l - 1], source_prod_matrix, np.array([0, 1, 0])
                    )

                if j == self.m and l == 0:
                    t1 = 0
                    t2 = self._FInt(
                        self.U[j - 1, l], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t3 = self._FInt(
                        self.L[j - 1, l], source_prod_matrix, np.array([0, 1, 0])
                    )
                    t4 = 0
                    t5 = 0
                    t6 = 0

                if j == self.m and l == self.n:
                    t1 = 0
                    t2 = 0
                    t3 = 0
                    t4 = self._FInt(
                        self.U[j - 1, l - 1], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t5 = 0
                    t6 = 0

                if j != 0 and j != self.m and l != 0 and l != self.n:
                    t1 = self._FInt(
                        self.L[j, l], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t2 = self._FInt(
                        self.U[j - 1, l], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t3 = self._FInt(
                        self.L[j - 1, l], source_prod_matrix, np.array([0, 1, 0])
                    )
                    t4 = self._FInt(
                        self.U[j - 1, l - 1], source_prod_matrix, np.array([1, 0, 0])
                    )
                    t5 = self._FInt(
                        self.L[j, l - 1], source_prod_matrix, np.array([0, 0, 1])
                    )
                    t6 = self._FInt(
                        self.U[j, l - 1], source_prod_matrix, np.array([0, 1, 0])
                    )

                # computing rh
                self.right_hand[idx] = t1 + t2 + t3 + t4 + t5 + t6

    def _FInt(self, T, fmatrix, v):
        # evaluating source term at f at the vertices of the element triangle
        f_11 = self._f_func(T["x"][0], T["y"][0], fmatrix)
        f_12 = self._f_func(T["x"][1], T["y"][1], fmatrix)
        f_13 = self._f_func(T["x"][2], T["y"][2], fmatrix)

        return self._trmatrix(T, f_11, f_12, f_13, v)

    def _f_func(self, x, y, f):
        """
        Evaluates the source term 'f' at each grid point
            % f is non-zero only at injection and production wells
            % x and y are the coordinates of the grid point
            % The corresponding index locations in the matrix for f
            % are determined in mm and nn respectively.
        """
        mm = int(round((x - self.left) / self.dx))
        nn = int(round((y - self.bottom) / self.dy))

        return f[nn, mm]

    def _trmatrix(self, T, f_1, f_2, f_3, v):
        s = self._polyarea(T["x"], T["y"])
        f_c = (f_1 + f_2 + f_3) / 3
        v_c = (v[0] + v[1] + v[2]) / 3
        f_4 = (f_2 + f_3) / 2
        f_5 = (f_1 + f_3) / 2
        f_6 = (f_1 + f_2) / 2

        v_4 = (v[1] + v[2]) / 2
        v_5 = (v[2] + v[0]) / 2
        v_6 = (v[0] + v[1]) / 2

        return (f_4 * v_4 + f_5 * v_5 + f_6 * v_6 + f_c * v_c) * s / 4

    def _set_A(self):
        """
        Assembles FEM stiffness matrix A
        Returns A as a sparse matrix

        Corresponds to the SetA.m function in matlab
        """
        self.A = np.zeros(((self.m + 1) * (self.n + 1) * 7, 3))
        list_idx = 0

        for j in range(self.m + 1):
            for l in range(self.n + 1):

                a = self.grid_size[j, l]
                idx = j + l * (self.m + 1)

                # center
                self.A[list_idx, :] = np.array([idx, idx, a["c"]])
                list_idx += 1

                # west
                if j > 0:
                    self.A[list_idx, :] = np.array([idx, idx - 1, a["w"]])
                    list_idx += 1

                # northwest
                if j > 0 and l < self.n:
                    self.A[list_idx, :] = np.array([idx, idx + self.m, a["nw"]])
                    list_idx += 1

                # north
                if l < self.n:
                    self.A[list_idx, :] = np.array([idx, idx + (self.m + 1), a["n"]])
                    list_idx += 1

                # east
                if j < self.m:
                    self.A[list_idx, :] = np.array([idx, idx + 1, a["e"]])
                    list_idx += 1

                # south
                if l > 0:
                    self.A[list_idx, :] = np.array([idx, idx - (self.m + 1), a["s"]])
                    list_idx += 1

                # southeast
                if j < self.m and l > 0:
                    self.A[list_idx, :] = np.array([idx, idx - self.m, a["se"]])
                    list_idx += 1

        self.A = self.A[:list_idx, :]

    def _set_B(self):
        self.B = np.zeros(((self.m + 1) * (self.n + 1), 1))

        for j in range(self.m + 1):
            for l in range(self.n + 1):

                a = self.grid_size[j, l]
                idx = j + l * (self.m + 1)

                self.B[idx] = a["const"]

        self.B = self.right_hand - self.B

    def get_A_B_matrices(self):
        self._set_A()
        self._set_B()

        # creating the sparsed matrix for A
        rows = np.ravel(self.A[:, 0].reshape(1, -1))
        cols = np.ravel(self.A[:, 1].reshape(1, -1))
        values = np.ravel(self.A[:, 2].reshape(1, -1))
        self.sparsed_A = coo_matrix(
            (values, (rows, cols)), shape=(np.shape(self.B)[0], np.shape(self.B)[0])
        )  # TODO: Verify correctedness of salf.sparsed_A
