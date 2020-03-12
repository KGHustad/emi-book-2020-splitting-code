import itertools

import numpy as np
import pytest

from index_mesh import IndexMesh

"""
@pytest.fixture
def m():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)
"""

def test_index_point():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)

    # x-dimension is contiguous and has stride 1, y-dimension has the middle stride, and the z-dimension has the longest stride
    count = 0
    for zi in range(Nz):
        for yi in range(Ny):
            for xi in range(Nx):
                assert m.index_point(xi, yi, zi) == count
                count += 1

def test_index_line_x():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)

    # x-dimension is contiguous and has stride 1, y-dimension has the middle stride, and the z-dimension has the longest stride
    for zi in range(Nz):
        for yi in range(Ny):
            offset = zi * Nx*Ny + yi * Nx

            # end-to-end
            #actual = m.index_line(0:Nx, yi, zi)
            actual = m[0:Nx, yi, zi]
            expected = np.arange(offset, offset + Nx)
            np.testing.assert_equal(actual, expected)

            middle = Nx//2

            # start-to-middle
            actual = m[0:middle, yi, zi]
            expected = np.arange(offset, offset + middle)
            np.testing.assert_equal(actual, expected)

            # middle-to-end
            actual = m[middle:Nx, yi, zi]
            expected = np.arange(offset + middle, offset + Nx)
            np.testing.assert_equal(actual, expected)

            # middle point
            actual = m[middle:middle+1, yi, zi]
            expected = np.arange(offset + middle, offset + middle+1)
            np.testing.assert_equal(actual, expected)

def test_index_line_y():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)

    # x-dimension is contiguous and has stride 1, y-dimension has the middle stride, and the z-dimension has the longest stride
    for zi in range(Nz):
        for xi in range(Nx):
            offset = zi * Nx*Ny + xi
            step = Nx

            # end-to-end
            actual = m[xi, 0:Ny, zi]
            expected = np.arange(offset, offset + Ny*step, step)
            np.testing.assert_equal(actual, expected)

            middle = Ny//2

            # start-to-middle
            actual = m[xi, 0:middle, zi]
            expected = np.arange(offset, offset + middle*step, step)
            np.testing.assert_equal(actual, expected)

            # middle-to-end
            actual = m[xi, middle:Ny, zi]
            expected = np.arange(offset + middle*step, offset + Ny*step, step)
            np.testing.assert_equal(actual, expected)

            # middle point
            actual = m[xi, middle:middle+1, zi]
            expected = np.arange(offset + middle*step, offset + (middle+1)*step, step)
            np.testing.assert_equal(actual, expected)

def test_index_line_z():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)

    # x-dimension is contiguous and has stride 1, y-dimension has the middle stride, and the z-dimension has the longest stride
    for yi in range(Nz):
        for xi in range(Nx):
            offset = yi * Nx + xi
            step = Nx*Ny

            # end-to-end
            actual = m[xi, yi, 0:Nz]
            expected = np.arange(offset, offset + Nz*step, step)
            np.testing.assert_equal(actual, expected)

            middle = Nz//2

            # start-to-middle
            actual = m[xi, yi, 0:middle]
            expected = np.arange(offset, offset + middle*step, step)
            np.testing.assert_equal(actual, expected)

            # middle-to-end
            actual = m[xi, yi, middle:Nz]
            expected = np.arange(offset + middle*step, offset + Nz*step, step)
            np.testing.assert_equal(actual, expected)

            # middle point
            actual = m[xi, yi, middle:middle+1]
            expected = np.arange(offset + middle*step, offset + (middle+1)*step, step)
            np.testing.assert_equal(actual, expected)

def test_index_plane_yz():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)
    mid_y = Ny//2
    mid_z = Nz//2
    for xi in range(Nx):
        # full plane
        actual = m[xi,:,:]
        expected = np.arange(xi, N, Nx)
        np.testing.assert_equal(actual, expected)

        # quarter plane (north, high)
        actual = m[xi,mid_y:,mid_z:]
        expected = np.zeros((Nz-mid_z, Ny-mid_y), dtype=np.int64)
        for i, zi in enumerate(range(mid_z, Nz)):
            for j, yi in enumerate(range(mid_y, Ny)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

        # quarter plane (south, low)
        actual = m[xi,:mid_y,:mid_z]
        expected = np.zeros((mid_z, mid_y), dtype=np.int64)
        for i, zi in enumerate(range(mid_z)):
            for j, yi in enumerate(range(mid_y)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

def test_index_plane_xz():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)
    mid_x = Nx//2
    mid_z = Nz//2
    for yi in range(Ny):
        # full plane
        start_x = 0
        stop_x = Nx
        start_z = 0
        stop_z = Nz

        actual = m[:,yi,:]
        expected = np.zeros((stop_z - start_z, stop_x - start_x), dtype=np.int64)
        for i, zi in enumerate(range(start_z, stop_z)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

        # quarter plane (east, high)
        start_x = mid_x
        stop_x = Nx
        start_z = mid_z
        stop_z = Nz

        actual = m[start_x:stop_x, yi, start_z:stop_z]
        expected = np.zeros((stop_z - start_z, stop_x - start_x), dtype=np.int64)
        for i, zi in enumerate(range(start_z, stop_z)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

        # quarter plane (west, low)
        start_x = 0
        stop_x = mid_x
        start_z = 0
        stop_z = mid_z

        actual = m[start_x:stop_x, yi, start_z:stop_z]
        expected = np.zeros((stop_z - start_z, stop_x - start_x), dtype=np.int64)
        for i, zi in enumerate(range(start_z, stop_z)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

def test_index_plane_xy():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)
    mid_x = Nx//2
    mid_y = Ny//2
    for zi in range(Nz):
        # full plane
        start_x = 0
        stop_x = Nx
        start_y = 0
        stop_y = Ny

        actual = m[:,:,zi]
        expected = np.zeros((stop_y - start_y, stop_x - start_x), dtype=np.int64)
        for i, yi in enumerate(range(start_y, stop_y)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

        # quarter plane (east, high)
        start_x = mid_x
        stop_x = Nx
        start_y = mid_y
        stop_y = Ny

        actual = m[start_x:stop_x, start_y:stop_y, zi]
        expected = np.zeros((stop_y - start_y, stop_x - start_x), dtype=np.int64)
        for i, yi in enumerate(range(start_y, stop_y)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

        # quarter plane (west, low)
        start_x = 0
        stop_x = mid_x
        start_y = 0
        stop_y = mid_y

        actual = m[start_x:stop_x, start_y:stop_y, zi]
        expected = np.zeros((stop_y - start_y, stop_x - start_x), dtype=np.int64)
        for i, yi in enumerate(range(start_y, stop_y)):
            for j, xi in enumerate(range(start_x, stop_x)):
                expected[i,j] = m.index_point(xi, yi, zi)
        expected = expected.flatten()
        np.testing.assert_equal(actual, expected)

def test_index_cuboid():
    Nx, Ny, Nz = 3, 4, 5
    N = Nx*Ny*Nz
    m = IndexMesh(Nx, Ny, Nz, N)
    mid_x = Nx//2
    mid_y = Ny//2
    mid_z = Nz//2

    # full
    actual = m[:, :, :]
    expected = np.arange(N, dtype=np.int64)
    np.testing.assert_equal(actual, expected)


    x_ranges = [(0, mid_x), (mid_x, Nx)]
    y_ranges = [(0, mid_y), (mid_y, Ny)]
    z_ranges = [(0, mid_z), (mid_z, Nz)]

    # test the eight cuboids that touch the middle and one of the mesh corners
    #for x_range, y_range, z_range in itertools.product([x_ranges, y_ranges, z_ranges]):
    for x_range, y_range, z_range in itertools.product(x_ranges, y_ranges, z_ranges):
        start_x, stop_x = x_range
        start_y, stop_y = y_range
        start_z, stop_z = z_range
        actual = m[start_x:stop_x, start_y:stop_y, start_z:stop_z]
        expected_3d = np.zeros((stop_z - start_z, stop_y - start_y, stop_x - start_x), dtype=np.int64)
        for i, zi in enumerate(range(start_z, stop_z)):
            for j, yi in enumerate(range(start_y, stop_y)):
                for k, xi in enumerate(range(start_x, stop_x)):
                    expected_3d[i,j,k] = m.index_point(xi, yi, zi)
        expected = expected_3d.flatten()
        np.testing.assert_equal(actual, expected)
