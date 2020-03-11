import math
import operator

import numpy as np
try:
    import numba
    import numba.typed
    has_numba = True
except ImportError:
    has_numba = False


if has_numba:
    @numba.njit(cache=True)
    def find_local_numbering_numba(N, local_numbering_array, input_arrays):
        """Finds the local numbering of each input array in local_numbering_array.
        The input arrays must be disjoint.
        """
        assert len(input_arrays) < 2**31-1
        vec = np.zeros(N, dtype=np.int32)

        # mark the input arrays in vec
        for i, input_array in enumerate(input_arrays):
            ind = i + 1
            vec[input_array] = ind

        output_arrays = [np.zeros_like(input_array) for input_array in input_arrays]
        output_array_positions = np.zeros(len(output_arrays), dtype=np.int64)
        # slide through vec
        for l, k in enumerate(local_numbering_array):
            val = vec[k]
            if val == 0:
                continue

            ind = val - 1
            pos = output_array_positions[ind]
            output_arrays[ind][pos] = l

            # advance position
            output_array_positions[ind] += 1

        for ind in range(len(output_arrays)):
            assert output_array_positions[ind] == output_arrays[ind].size

        return output_arrays

def find_local_numbering_python(N, local_numbering_array, input_arrays):
    """Finds the local numbering of each input array in local_numbering_array.
    The input arrays must be disjoint.
    """
    assert len(input_arrays) < 2**31-1
    vec = np.zeros(N, dtype=np.int32)

    # mark the input arrays in vec
    for i, input_array in enumerate(input_arrays):
        ind = i + 1
        vec[input_array] = ind

    output_arrays = [np.zeros_like(input_array) for input_array in input_arrays]
    output_array_positions = np.zeros(len(output_arrays), dtype=np.int64)
    # slide through vec
    for l, k in enumerate(local_numbering_array):
        val = vec[k]
        if val == 0:
            continue

        ind = val - 1
        pos = output_array_positions[ind]
        output_arrays[ind][pos] = l

        # advance position
        output_array_positions[ind] += 1

    for ind in range(len(output_arrays)):
        assert output_array_positions[ind] == output_arrays[ind].size

    return output_arrays

def find_local_numbering(N, local_numbering_array, input_arrays):
    if has_numba:
        input_arrays_typed = numba.typed.List()
        [input_arrays_typed.append(a) for a in input_arrays]
        return find_local_numbering_numba(N, local_numbering_array, input_arrays_typed)
    else:
        return find_local_numbering_python(N, local_numbering_array, input_arrays)

def select_best_opencl_device():
    import pyopencl
    platforms = pyopencl.get_platforms()

    # find all devices that support FP64
    fp64_devices = []
    for platform in platforms:
        devices = platform.get_devices()
        for device in devices:
            if device.native_vector_width_double != 0:
                fp64_devices.append(device)

    if len(fp64_devices) == 0:
        raise ValueError("Could not find any FP64-capable devices")

    # prefer GPU devices
    gpu_devices = list(filter(
        lambda device: device.type == pyopencl.device_type.GPU,
        fp64_devices
    ))
    if len(gpu_devices) == 0:
        # just go with device 0, which is usually a CPU
        return fp64_devices[0]
    elif len(gpu_devices) == 1:
        # choose the only FP64-capable GPU
        return gpu_devices[0]

    # in lack of better heuristics, pick the GPU with the most compute units
    ranked_gpu_devices = sorted(gpu_devices, key=operator.attrgetter('max_compute_units'), reverse=True)
    return ranked_gpu_devices[0]

def ceil_to_multiple(a, d):
    return (a+d-1) // d * d

def round_to_nearest(x):
    return math.floor(x+0.5)
