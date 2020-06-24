import numpy as np
import ComPASS


def compare_dumps(a):
    print("From Python:")
    print(a)
    print("From Fortran:")
    ComPASS.dump_array_in_fortran(a)


def test_dump_array():
    # FIXME: this is due to the fact each eos embeds everything
    ComPASS.load_eos("linear_water")
    a = np.ascontiguousarray(np.arange(6), dtype=np.double)
    a.shape = (-1, 2)
    compare_dumps(a)
    ComPASS.increment_first_column_in_fortran(a)
    compare_dumps(a)


def test_increment_array():
    # FIXME: this is due to the fact each eos embeds everything
    ComPASS.load_eos("linear_water")
    a = np.ascontiguousarray(np.random.random((3, 2)), dtype=np.double)
    b = np.copy(a)
    ComPASS.increment_first_column_in_fortran(a)
    b[:, 0] += 1
    assert np.allclose(a, b)


if __name__ == "__main__":
    test_increment_array()
