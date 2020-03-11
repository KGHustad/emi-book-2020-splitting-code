import numpy as np

empty_array = np.zeros(0, dtype=np.int64)

class IndexMesh(object):
    def __init__(self, Nx, Ny, Nz, N):
        assert Nx*Ny*Nz == N
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.N = N

    def __getitem__(self, idx):
        assert len(idx) == 3 # mesh is 3D
        xi, yi, zi = idx
        arg_is_slice = [type(arg) == slice for arg in (xi, yi, zi)]
        num_slices = sum(arg_is_slice)
        if num_slices == 0:
            return self.index_point(xi, yi, zi)
        elif num_slices == 1:
            return self.index_line(xi, yi, zi)
        elif num_slices == 2:
            return self.index_plane(xi, yi, zi)
        elif num_slices == 3:
            return self.index_cuboid(xi, yi, zi)
        else:
            raise ValueError(f"Invalid number of slices: {num_slices}")


    def index_point(self, xi, yi, zi):
        # raise an error if xi, yi or zi are out of bounds
        # don't raise an error if the values are one larger than the maximum valid, since this is often used to specify end points
        if xi < 0 or self.Nx < xi:
            raise ValueError("xi = %d exceeds the mesh interval [%d, %d]" % (xi, 0, self.Nx))
        if yi < 0 or self.Ny < yi:
            raise ValueError("yi = %d exceeds the mesh interval [%d, %d]" % (yi, 0, self.Ny))
        if zi < 0 or self.Nz < zi:
            raise ValueError("zi = %d exceeds the mesh interval [%d, %d]" % (zi, 0, self.Nz))

        # x-dimension is contiguous and has stride 1, y-dimension has the middle stride, and the z-dimension has the longest stride
        return xi + yi*self.Nx + zi*self.Nx*self.Ny

    def index_line(self, xi, yi, zi):
        Nx = self.Nx
        Ny = self.Ny
        Nz = self.Nz

        arg_is_slice = [type(arg) == slice for arg in (xi, yi, zi)]
        assert sum(arg_is_slice) == 1
        if type(xi) == slice:
            #offset = self.index_point(0, yi, zi)
            start, stop, step = xi.indices(Nx)
            assert step == 1

            #return np.arange(start+offset, stop+offset)
            start_1d = self.index_point(start, yi, zi)
            stop_1d = self.index_point(stop, yi, zi)
            step_1d = self.index_point(start+1, yi, zi) - start_1d

            return np.arange(start_1d, stop_1d, step_1d)
        if type(yi) == slice:
            start, stop, step = yi.indices(Ny)
            assert step == 1

            start_1d = self.index_point(xi, start, zi)
            stop_1d = self.index_point(xi, stop, zi)
            step_1d = self.index_point(xi, start+1, zi) - start_1d

            return np.arange(start_1d, stop_1d, step_1d)
        if type(zi) == slice:
            start, stop, step = zi.indices(Nz)
            assert step == 1

            start_1d = self.index_point(xi, yi, start)
            stop_1d = self.index_point(xi, yi, stop)
            step_1d = self.index_point(xi, yi, start+1) - start_1d

            return np.arange(start_1d, stop_1d, step_1d)
        raise ValueError("Found no slice. This should never occur.")

    def index_plane(self, xi, yi, zi):
        Nx = self.Nx
        Ny = self.Ny
        Nz = self.Nz

        arg_is_slice = [type(arg) == slice for arg in (xi, yi, zi)]
        assert sum(arg_is_slice) == 2
        if type(xi) is not slice:
            # plane in yz

            start_y, stop_y, step_y = yi.indices(Ny)
            start_z, stop_z, step_z = zi.indices(Nz)
            assert step_y == 1
            assert step_z == 1
            empty_range = not (start_y < stop_y and start_z < stop_z)
            if empty_range:
                return empty_array
            assert start_y < stop_y
            assert start_z < stop_z


            start_1d = self.index_point(xi, start_y, start_z)
            #stop_1d = self.index_point(xi, stop_y, stop_z)
            step_1d_y = self.index_point(xi, start_y+1, start_z) - start_1d
            step_1d_z = self.index_point(xi, start_y, start_z+1) - start_1d

            y_len = max((stop_y - start_y - 1) // step_y, 0) + 1
            z_len = max((stop_z - start_z - 1) // step_z, 0) + 1
            if step_1d_y < step_1d_z:
                indices_2d = np.zeros((z_len, y_len), dtype=np.int64)
                index_line_z = self.index_line(xi, 0, zi)
                index_line_y = self.index_line(0, yi, 0)

                indices_2d += index_line_z[:,np.newaxis]
                indices_2d += index_line_y[:,np.newaxis].T

                return indices_2d.flatten()
            else:
                indices_2d = np.zeros((y_len, z_len), dtype=np.int64)
                index_line_y = self.index_line(xi, yi, 0)
                index_line_z = self.index_line(0, 0, zi)

                indices_2d += index_line_y[:,np.newaxis]
                indices_2d += index_line_z[:,np.newaxis].T

                return indices_2d.flatten()

        if type(yi) is not slice:
            # plane in xz

            start_x, stop_x, step_x = xi.indices(Nx)
            start_z, stop_z, step_z = zi.indices(Nz)
            assert step_x == 1
            assert step_z == 1
            empty_range = not (start_x < stop_x and start_z < stop_z)
            if empty_range:
                return empty_array
            assert start_x < stop_x
            assert start_z < stop_z

            start_1d = self.index_point(start_x, yi, start_z)
            #stop_1d = self.index_point(stop_x, yi, stop_z)
            step_1d_x = self.index_point(start_x+1, yi, start_z) - start_1d
            step_1d_z = self.index_point(start_x, yi, start_z+1) - start_1d

            x_len = max((stop_x - start_x - 1) // step_x, 0) + 1
            z_len = max((stop_z - start_z - 1) // step_z, 0) + 1
            if step_1d_x < step_1d_z:
                indices_2d = np.zeros((z_len, x_len), dtype=np.int64)
                index_line_z = self.index_line(0, yi, zi)
                index_line_x = self.index_line(xi, 0, 0)

                indices_2d += index_line_z[:,np.newaxis]
                indices_2d += index_line_x[:,np.newaxis].T

                return indices_2d.flatten()
            else:
                indices_2d = np.zeros((x_len, z_len), dtype=np.int64)
                index_line_x = self.index_line(xi, yi, 0)
                index_line_z = self.index_line(0, 0, zi)

                indices_2d += index_line_x[:,np.newaxis]
                indices_2d += index_line_z[:,np.newaxis].T

                return indices_2d.flatten()
        if type(zi) is not slice:
            # plane in xy

            start_x, stop_x, step_x = xi.indices(Nx)
            start_y, stop_y, step_y = yi.indices(Ny)
            assert step_x == 1
            assert step_y == 1
            empty_range = not (start_x < stop_x and start_y < stop_y)
            if empty_range:
                return empty_array
            assert start_x < stop_x
            assert start_y < stop_y

            start_1d = self.index_point(start_x, start_y, zi)
            #stop_1d = self.index_point(stop_x, stop_y, zi)
            step_1d_x = self.index_point(start_x+1, start_y, zi) - start_1d
            step_1d_y = self.index_point(start_x, start_y+1, zi) - start_1d

            x_len = max((stop_x - start_x - 1) // step_x, 0) + 1
            y_len = max((stop_y - start_y - 1) // step_y, 0) + 1
            if step_1d_x < step_1d_y:
                indices_2d = np.zeros((y_len, x_len), dtype=np.int64)
                index_line_y = self.index_line(0, yi, zi)
                index_line_x = self.index_line(xi, 0, 0)

                indices_2d += index_line_y[:,np.newaxis]
                indices_2d += index_line_x[:,np.newaxis].T

                return indices_2d.flatten()
            else:
                indices_2d = np.zeros((x_len, y_len), dtype=np.int64)
                index_line_x = self.index_line(xi, 0, zi)
                index_line_y = self.index_line(0, yi, 0)

                indices_2d += index_line_x[:,np.newaxis]
                indices_2d += index_line_y[:,np.newaxis].T

                return indices_2d.flatten()
        raise ValueError("Did not find the right plane. This should never occur.")

    def index_cuboid(self, xi, yi, zi):
        arg_is_slice = [type(arg) == slice for arg in (xi, yi, zi)]
        assert sum(arg_is_slice) == 3

        Nx = self.Nx
        Ny = self.Ny
        Nz = self.Nz

        start_x, stop_x, step_x = xi.indices(Nx)
        start_y, stop_y, step_y = yi.indices(Ny)
        start_z, stop_z, step_z = zi.indices(Nz)
        assert step_x == 1
        assert step_y == 1
        assert step_z == 1

        x_len = max((stop_x - start_x - 1) // step_x, 0) + 1
        y_len = max((stop_y - start_y - 1) // step_y, 0) + 1
        z_len = max((stop_z - start_z - 1) // step_z, 0) + 1

        indices_3d = np.zeros((z_len, y_len, x_len), dtype=np.int64)
        index_line_x = self.index_line(xi, 0, 0)
        index_line_y = self.index_line(0, yi, 0)
        index_line_z = self.index_line(0, 0, zi)

        indices_3d += index_line_x[np.newaxis,np.newaxis,:]
        indices_3d += index_line_y[np.newaxis,:,np.newaxis]
        indices_3d += index_line_z[:,np.newaxis,np.newaxis]
        return indices_3d.flatten()

