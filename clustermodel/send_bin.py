import numpy as np

def send_array(arr, fh):
    # number of probes (columns)
    arr = np.asarray(arr).T
    shape = arr.shape
    if len(shape) == 1:  # could be e.g.: (24,). we want (24, 1)
        shape = (shape[0], 1)

    np.array(shape, dtype=np.int64).tofile(fh)
    np.asarray(arr).flatten().astype(np.float64).tofile(fh)

def send_arrays(arrays, fh):
    """
    format is to send:
        + int64 of number of matrices
        + for each matrix:
            + int64, int64 of shape
            + float64 * (int64 * int64) of values
    """
    fh.seek(0)
    np.array([len(arrays)], dtype=np.int64).tofile(fh)
    for array in arrays:
        send_array(array, fh)
    fh.flush()


if __name__ == "__main__":

    fh = open('t.bin', 'w')
    a1 = np.arange(100).reshape((20, 5)) + 0.11
    a2 = np.arange(46).reshape((23, 2)) + 0.22
    print a2
    a3 = np.arange(88).reshape((8, 11)) + 0.33
    send_arrays((a1, a2, a3), fh)


