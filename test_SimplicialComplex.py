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
    L1 = octahedron.link([1])
    assert L1.f_vector() == (1,4,4)
    assert set(L1.facets()) == set([(2,3), (3,5), (5,6), (2,6)])

    L2 = non_pure_complex.link([3])
    assert set(L2.facets()) == set([(1,2), (2,4), (5,)])
    assert not L2.is_pure()

def test_add_face():
    K = SimplicialComplex([[1,2,3], [1,3,5], [1,5,6], [1,2,6],
                                    [2,3,4], [3,4,5], [4,5,6], [2,4,6]])
    K.add_face([1,4])
    assert K.f_vector() == (1,6,13,8)

def test_remove_face():
    K = SimplicialComplex([[1,2,3], [1,3,5], [1,5,6], [1,2,6],
                                    [2,3,4], [3,4,5], [4,5,6], [2,4,6]])
    K.add_face([1,4])
    K.remove_face([1,4])
    assert set(K.facets()) == set(octahedron.facets())

    K.add_face([1,2,3,4])
    #assert set(K.facets()) == set([[1,2,3,4], [1,3,5], [1,5,6], [1,2,6], [3,4,5], [4,5,6], [2,4,6]])


def test_deletion():
    D1 = octahedron.deletion([1])
    assert octahedron.f_vector() == (1, 6, 12, 8)
    assert set(D1.facets()) == {(2,3,4), (3,4,5), (4,5,6), (2,4,6)}

    D1 = D1.deletion([4])
    # This causes the dimension to drop
    assert set(D1.facets()) == {(2,3), (3,5), (5,6), (2,6)}
    assert D1.get_dimension() == 1
    assert D1.get_rank() == 2
    assert set(D1.vertices()) == {2,3,5,6}

def test_remove_face():
    # K is a copy of the octahedron
    K = SimplicialComplex([[1,2,3], [1,3,5], [1,5,6], [1,2,6],
                                    [2,3,4], [3,4,5], [4,5,6], [2,4,6]])

    # Test that removing a non-existent face does nothing.
    K.remove_face([1,3,7])
    assert K.f_vector() == (1,6,12,8)

    # Removing vertex 1 results in a cone over a 4-cycle with apex vertex 4.
    # It has 4 facets and 5 vertices.  Vertex 1 has been removed from the
    # vertex set.
    K.remove_face([1])
    assert set(K.facets()) == {(2,3,4), (3,4,5), (4,5,6), (2,4,6)}
    assert K.vertices() == {2,3,4,5,6}

    # Removing vertex 4 results in a 4-cycle. The dimension has dropped from
    # 2 to 1. There are now only four vertices.
    K.remove_face([4])
    assert K.get_dimension() == 1
    assert K.get_rank() == 2
    assert K.vertices() == {2, 3, 5, 6}

    # Removing edge {2,3} gives a path on four vertices.  The complex is still
    # one-dimensional, and it still has four vertices.
    K.remove_face([2,3])
    assert K.get_dimension() == 1
    assert K.get_rank() == 2
    assert K.f_vector() == (1, 4, 3)
    assert K.vertices() == {2, 3, 5, 6}

    # Now we will remove three vertices, 2, 3, and 5. This leaves a complex
    # consisting of a single vertex.  The complex is now zero-dimensional.
    K.remove_face([2])
    K.remove_face([3])
    K.remove_face([5])
    assert K.get_dimension() == 0
    assert K.f_vector() == (1,1)

    # Finally, removing vertex 6 results in the complex whose only face is the
    # empty set. The dimension is -1 and its vertex set is empty.
    K.remove_face([6])
    assert K.get_dimension() == -1
    assert K.vertices() == set()

    # What happens if we try to remove the empty face?
    L = SimplicialComplex([[1,2,3],[2,3,4]])
    L.remove_face([])
    assert L.get_rank() == 0



if __name__ == '__main__':
    K = SimplicialComplex([[1,2,3], [1,3,5], [1,5,6], [1,2,6],
                                    [2,3,4], [3,4,5], [4,5,6], [2,4,6]])
    print(K.deletion([1]))
    K.add_face([1,2,3,4])
    print(K.f_vector())
    print(K)
