"""
Testing suite for SimplicialComplex package
"""

from SimplicialComplex import *



octahedron = SimplicialComplex([[1,2,3], [1,3,5], [1,5,6], [1,2,6],
                                [2,3,4], [3,4,5], [4,5,6], [2,4,6]])

non_pure_complex = SimplicialComplex([(1,2,3), (2,3,4), (3,5), (4,5)])

mobius = SimplicialComplex([[1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6],
                            [5, 6, 1], [6, 1, 2]])
projective_plane = SimplicialComplex([
                [1, 2, 3], [1, 4, 5], [1, 5, 6], [1, 2, 6], [2, 4, 6],
                [2, 4, 5], [2, 3, 5], [3, 5, 6], [3, 4, 6], [1, 3, 4]])

torus = SimplicialComplex([
    [1,2,5], [2,4,5], [1,4,6], [1,2,6], [4,6,7], [2,4,7], [2,3,7],
    [1,3,7], [1,5,7], [5,6,7], [3,5,6], [3,4,5], [1,3,4], [2,3,6]
])

empty_complex = SimplicialComplex()

def test_is_pure():
    assert octahedron.is_pure()
    assert not non_pure_complex.is_pure()
    assert torus.is_pure()
    assert projective_plane.is_pure()
    assert mobius.is_pure()
    assert empty_complex.is_pure()

def test_f_vector():
    assert octahedron.f_vector() == (1, 6, 12, 8)
    assert torus.f_vector() == (1, 7, 21, 14)
    assert empty_complex.f_vector() == (1,)
    assert non_pure_complex.f_vector() == (1, 5, 7, 2)

def test_h_vector():
    assert octahedron.h_vector() == (1, 3, 3, 1)
    assert torus.h_vector() == (1, 4, 10, -1)
    assert empty_complex.h_vector() == (1, )
    assert non_pure_complex.h_vector() == (1, 2, 0, -1)

def test_get_dimension():
    assert octahedron.get_dimension() == 2
    assert torus.get_dimension() == 2
    assert empty_complex.get_dimension() == -1
    assert non_pure_complex.get_dimension() == 2

def test_get_rank():
    assert octahedron.get_rank() == 3
    assert torus.get_rank() == 3
    assert empty_complex.get_rank() == 0
    assert non_pure_complex.get_rank() == 3

def test_euler_characteristic():
    assert octahedron.euler_characteristic() == 1
    assert octahedron.euler_characteristic(reduced = False) == 2
    assert octahedron.euler_characteristic(reduced = True) == 1

def test_link():
    L = octahedron.link([1])
    assert L.f_vector() == (1,4,4)
