import numpy as np
from scipy.sparse import coo_matrix
from matplotlib.tri import Triangulation


class Grid:
    """
    Encapsulates mesh generation, triangulation, FEM matrix assembly, and source vector setup.
    """
    def __init__(self, m: int = 32, n: int = 32, left: float = 0.0, right: float = 1.0, bottom: float = 0.0, top: float = 1.0):
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
    
    ## FIXME: Functions below need to be updated
    def set_FE_meshgrid(self):
        """
        Generate FE coordinate grid for elliptic pressure calculations
        Analogous to the setGrid.m function in the MATLAB code
        """
        pass

    def set_triangulation(self):
        xv = self.x.flatten()
        yv = self.y.flatten()
        self.tri = Triangulation(xv, yv)

    def set_right_hand(self, rhs_func):
        """
        Sets the right-hand side (source) from a function.
        Args:
            rhs_func: takes (x, y) arrays and returns a same-shaped source field
        """
        self.RHS = rhs_func(self.x, self.y)

    def set_A(self, beta_field: np.ndarray):
        """
        Assembles FEM stiffness matrix A
        """
        rows, cols, data = [], [], []
        num_nodes = (self.m + 1) * (self.n + 1)

        def idx(i, j): return i * (self.m + 1) + j

        for i in range(self.n + 1):
            for j in range(self.m + 1):
                center = idx(i, j)

                if i > 0:
                    up = idx(i - 1, j)
                    rows.append(center)
                    cols.append(up)
                    data.append(-beta_field[i, j] / self.dy ** 2)
                if i < self.n:
                    down = idx(i + 1, j)
                    rows.append(center)
                    cols.append(down)
                    data.append(-beta_field[i, j] / self.dy ** 2)
                if j > 0:
                    left = idx(i, j - 1)
                    rows.append(center)
                    cols.append(left)
                    data.append(-beta_field[i, j] / self.dx ** 2)
                if j < self.m:
                    right = idx(i, j + 1)
                    rows.append(center)
                    cols.append(right)
                    data.append(-beta_field[i, j] / self.dx ** 2)

                rows.append(center)
                cols.append(center)
                data.append(2 * beta_field[i, j] * (1 / self.dx ** 2 + 1 / self.dy ** 2))

        self.A = coo_matrix((data, (rows, cols)), shape=(num_nodes, num_nodes)).tocsc()

    def set_B(self, source_array: np.ndarray | None=None):
        """
        Sets vector B (Right Hand Side of linear system).
        """
        num_nodes = (self.m + 1) * (self.n + 1)
        self.B = np.zeros(num_nodes)
        if source_array is not None:
            self.B[:] = source_array.flatten()

    def get_flat_index_matrix(self) -> np.ndarray:
        """
        Returns a matrix of shape (n+1, m+1) with flat indices at each grid point.
        """
        return np.arange((self.m + 1) * (self.n + 1)).reshape((self.n + 1, self.m + 1))

