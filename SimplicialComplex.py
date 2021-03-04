"""
File: SimplicialComplex.py
Class: SimplicialComplex

Defines a simplicial complex object.  A simplicial complex is defined in terms
of two pieces of information:

1. A set, V, of vertices
2. A set of /faces/: a collection of subsets of V that is
   /closed under inclusion/, meaning that if F is a face and G is a subset of F,
   then G is a face as well.
"""

class Simplex:
    """
    Defines a Simplex object to be used in constructing a simplicial complex

    A simplex is a face in a simplicial complex. It consists of a set of
    vertices (elements), together with all possible subsets of that set of
    vertices as faces.

    Attributes
    ----------
    face_dictionary: dictionary
        A dictionary storing the faces of a simplex

    elements: tuple
        Tuple of vertices that form the simplex.

    rank: int
        The cardinality of the simplex.

    dimension: int
        The dimension of the simplex; one less than its rank.

    Methods
    -------
    faces: dictionary
        Returns the face dictionary.

    get_dimension: int
        Returns the dimension of the simplex.

    get_rank: int
        Returns the rank (cardinality) of the simplex.

    """

    def __init__(self, elements):
        """
        Parameters
        ----------
        elements: tuple, list, or set
            The vertices that comprise the simplex
            Will be converted to a tuple upon instantiation
        """

        # Store elements as a tuple
        self.elements = tuple(elements)

        # rank = size of the simplex
        self.rank = len(self.elements)
        self.dimension = self.rank - 1

        # Store the faces of a simplex as a dictionary where
        # keys = dimensions (-1 through self.dimension)
        # values = set of faces, stored as tuples
        self.face_dictionary = {}

        # Initially, the only face is the simplex itself.
        self.face_dictionary[self.dimension] = {self.elements}

        # Work through possible dimensions of faces, from self.dimension
        # down to -1.  Faces of dimension k are obtained from faces of dimension
        # k + 1 by removing one vertex at a time.
        # current_dim is used to track the dimension of faces that are currently
        # being added to the dictionary
        current_dim = self.dimension-1
        while current_dim >= -1:
            # Start with an empty set of faces of dimension current_dim
            self.face_dictionary[current_dim] = set([])
            # For each face of dimension one larger than current_dim and
            # each vertex v in that face, remove v to get a face of dimension
            # current_dimension.  Add it to the dictionary if it is not already
            # present.
            for existing_face in self.face_dictionary[current_dim+1]:
                for v in existing_face:
                    new_face = tuple([u for u in existing_face if u != v])
                    self.face_dictionary[current_dim].add(new_face)
            # Decrement the dimension by 1 after all faces have been tested.
            current_dim -= 1

    def __str__(self):
        return str(f'Simplex on vertex set {self.elements}')


    def faces(self):
        """Returns the dictionary of faces in the simplex."""
        return self.face_dictionary

    def get_dimension(self):
        """Returns the dimension of the simplex (size minus 1)."""
        return int(self.dimension)

    def get_rank(self):
        """Returns the size of the simplex."""
        return self.rank


###########################

def binomial(n,k):
    """ Returns the value of the binomial coefficient n choose k

    The binomial coefficient n choose k counts the number of k-element subsets
    of an n-element set.

    Parameters
    ----------
    n: int
        The size of a universal set of objects

    k: int
        The size of subsets to select from the universal set.
    """

    # If k < 0 or k > n, then n choose k = 0
    if k < 0 or k > n:
        return 0

    # If k == 0 or k == n, then n choose k = 1
    elif k == 0 or k == n:
        return 1

    # Otherwise, we can recursively compute n choose k using Pascal's relation
    # n choose k = (n-1 choose k-1) + (n-1 choose k)
    else:
        return binomial(n-1, k-1) + binomial(n-1, k)

###########################

class SimplicialComplex:
    """Defines a simplicial complex object.

    A simplicial complex is a collection of sets that is closed under inclusion.
    These sets are called /faces/. Each face is a Simplex.

    Attributes
    ----------
    face_dictionary: dictionary
        A dictionary storing the faces of the simplicial complex.

    dimension: int
        The dimension of the simplicial complex (max dimension of its faces).

    rank: int
        The rank of the simplicial complex (max rank of its faces).

    vertex_set: set
        The set of vertices of the simplicial complex.  We use the convention
        that the vertices are elements that appear in a face of the complex.
        Loops are not allowed.

    Methods
    -------
    get_dimension: int
        Returns the dimension of the complex

    get_rank: int
        Returns the rank of the complex

    faces: dictionary
        Returns the face dictionary

    vertices: set
        Returns the vertex set

    f-vector: tuple
        Returns the f-vector of the complex.


    h-vector: tuple
        Returns the h-vector of the complex.


    euler_characteristic: int, optional
        Returns the Euler characteristic of the complex. Default is to return
        the /reduced/ Euler characteristic.

    facets: set
        Returns the set of facets in the complex.

    is_pure: boolean
        Tests whether the complex is pure. A complex is /pure/ if all of its
        facets have the same dimension.

    add_face(new_face): self
        Adds the face new_face to self.

    link(F): SimplicialComplex
        Return the link of face F in the given simplicial complex.

    remove_face(F): self
        Removes face F from the given simplicial complex.

    deletion(F): SimplicialComplex
        Return the simplicial complex obtained by deleting face F from the
        simplicial complex. This does not change the simplicial complex itself.

    """

    def __init__(self, generating_faces = []):
        """
        A SimplicialComplex object is defined from a list of given faces, which
        are generally assumed to be facets (maximal faces under inclusion). A
        simplicial complex is generated from a list of faces by taking all
        possible subsets of those faces so that the resulting set family is
        closed under inclusion.

        If generating_faces is empty, the resulting simplicial complex consists
        only of the empty face.  We do not allow for the empty complex, which
        has no faces.

        Parameters
        ----------

        generating_faces: list (default = [])
            A list of faces used to define the complex.
        """

        # The rank of a complex is the size of the largest facet in the list of
        # generators. If generating_faces is empty, the resulting complex has
        # rank 0 by default.
        try:
            self.rank = max([len(S) for S in generating_faces])
        except:
            self.rank = 0

        # The dimension of a simplicial complex is one less than its rank.
        self.dimension = self.rank - 1

        # Simultaneously determine the vertex set and face dictionary of the
        # complex. In the face dictionary, keys are dimensions of faces, and
        # values are sets of tuples, the faces of each given dimension.
        # Initially, for each possible dimension of face, the set of faces of
        # that dimension is empty.
        self.vertex_set = set([])
        self.face_dictionary = {k: set() for k in range(-1,self.rank)}

        # The only face of dimension -1 is the empty face.
        self.face_dictionary[-1] = {tuple()}

        # For each face in the list of facets, create the simplex on its set of
        # vertices, and add all of those faces to the face dictionary. Add its
        # vertices to the vertex set.
        for new_facet in generating_faces:
            self.vertex_set = self.vertex_set.union(set(new_facet))
            # Get face dictionary from Simplex object generated by the vertices
            # in new_facet. Then add those faces to the face dictionary of the
            # simplicial complex. Faces are sorted to prevent adding multiple
            # instances of the same face.
            new_facet = sorted(new_facet)
            new_simplex = Simplex(new_facet).faces()
            for dimension in new_simplex.keys():
                for simplex_face in new_simplex[dimension]:
                    self.face_dictionary[dimension].add(simplex_face)


    def __str__(self):
        output_str = f'Simplicial complex of dimension {self.dimension} ' \
                     + f'with {len(self.vertex_set)} vertices and ' \
                     + f'{len(self.facets())} facets'
        return output_str

    def __repr__(self):
        return str(self)

    def get_dimension(self):
        """ Returns the dimension of the simplicial complex.

        The dimension of a face in a simplicial complex is its size minus one.
        The dimension of a simplicial complex is the max dimension of any of
        its faces.
        """

        return self.dimension

    def get_rank(self):
        """ Returns the rank of the simplicial complex.

        The rank of a face in a simplicial complex is its size. The rank of a
        simplicial complex is the max rank of any of its faces.
        """

        return self.rank

    def faces(self):
        """Returns the dictionary of faces in the simplicial complex.

        The faces are returned as a dictionary where keys are dimensions of
        faces and values are a set of tuples of faces of that dimension.
        """

        return self.face_dictionary

    def vertices(self):
        """Returns the vertex set of the simplicial complex.

        By construction, only vertices that appear in some face are allowed.
        In other words, loops are not allowed.
        """

        return self.vertex_set

    def f_vector(self):
        """ Returns the f-vector of the simplicial complex.

        The f-vector is the vector
            (f_{-1}, f_0, ..., f_{d-1}),
        where d-1 is the dimension of the complex and f_i is the number of
        i-dimensional faces in the complex.
        Note: Indexing in the f-vector is shifted by one, so f_vec[i] is the
        number of faces of size i, which is the number of (i-1)-dimensional
        faces.  For example, f_vec[1] is the number of vertices (0-dimensional
        faces).
        """

        f_vec = [len(self.face_dictionary[k]) for k in range(-1,self.rank)]
        return tuple(f_vec)

    def h_vector(self):
        """Returns the h-vector of the simplicial complex.

        The h-vector is the vector
            (h_0, h_1, ..., h_d),
        where d-1 is the dimension of the complex and h_j is defined by
            h_j = sum((-1)**(j-i) * binom(d-i, d-j) * f_{i-1} for 0 <= i <= j)
        """

        f_vec = self.f_vector()
        h_vec = tuple()
        # Append the h-numbers h_j one at a time for 0 <= j <= d.
        for j in range(0,self.rank+1):
            h_j = 0
            for i in range(0,j+1):
                h_j += (-1)**(j-i) * binomial(self.rank-i, self.rank-j) * f_vec[i]
            h_vec = (*h_vec, h_j)
        return h_vec

    def euler_characteristic(self, reduced = True):
        """Returns the Euler characteristic of a simplicial complex. By default,
        the reduced Euler characteristic is returned.

        Euler characteristic is defined as the alternating sum of the number
        of faces in each dimension:
            -f_{-1} + f_0 - f_1 + f_2 - f_3 + ... + (-1)^(d-1) * f_{d-1}
        The reduced Euler characteristic includes the empty face.  Non-reduced
        Euler characteristic does not.
        """
        euler_char = 0
        f_vec = self.f_vector()
        for i in range(self.rank + 1):
            euler_char += ((-1)**(i-1)) * f_vec[i]
        if reduced:
            return int(euler_char)
        else:
            # For non-reduced Euler characteristic, we need to exclude the -1
            # coming from (-1)**(-1) * f_{-1} from the existing sum. We do this
            # by simply adding one.
            return int(euler_char + 1)

    def facets(self):
        """Returns the set of facets of a simplicial complex.

        Facets are maximal faces under inclusion. They may be different than
        the set of faces used to define the complex.
        """

        # Every face of maximal dimension is a facet
        facet_list = list(self.face_dictionary[self.dimension])

        # Work through the facet dictionary, starting with codimension-one
        # faces, until we reach the empty face. Any face not contained in an
        # existing facet is a new facet.
        current_dimension = self.dimension - 1
        while current_dimension > -1:
            for test_face in self.face_dictionary[current_dimension]:
                # If the face is not contained in any existing facets,
                # add it to the facet list.
                is_contained = any(set(test_face).issubset(set(F))
                                   for F in facet_list)
                if not(is_contained):
                    facet_list.append(test_face)
            current_dimension -= 1
        return facet_list

    def is_pure(self):
        """Checks to see if the simplicial complex is pure.

        A simplicial complex is pure if all of its facets have the same
        dimension.
        """

        # The complex is not pure if there exists a facet whose cardinality
        # is less than the rank of the complex. Otherwise, it is pure.
        for facet in self.facets():
            if len(facet) != self.rank:
                return False
        return True


    def add_face(self, new_face):
        """Adds the given new face to the simplicial complex.
        """
        new_face = tuple(sorted(new_face))
        new_simplex = Simplex(new_face)
        for i in range(-1, new_simplex.get_rank()):
            for face in new_simplex.faces()[i]:
                # It could happen that the new face has larger dimension than
                # the given complex. In this case, we need to add new keys to
                # the dictionary before adding the new faces.
                if i not in self.face_dictionary.keys():
                    self.face_dictionary[i] = set([face])
                else:
                    self.face_dictionary[i].add(face)

        # Update the dimension, rank, and vertex_set attributes of the
        # simplicial complex to account for the fact that we may have
        # changed the dimension or increased the number of vertices.
        self.dimension = max(self.face_dictionary.keys())
        self.rank = self.dimension + 1
        new_vertices = [v for v in new_face if v not in self.vertex_set]
        self.vertex_set = self.vertex_set.union(new_vertices)

    def link(self, F):
        """Returns the link of face F in the given simplicial complex.

        The link of face F is defined to be the set of all faces in the given
        simplicial complex that contain F. It is a simplicial complex on its
        own.
        """
        # The link is only defined if F is a face of the simplicial complex.
        F = list(F)
        F = sorted(F)
        F = tuple(F)
        if F not in self.face_dictionary[len(F) - 1]:
            print(f'Face {F} does not belong to the simplicial complex.')
            return None

        else:
            # The facets of link(F) are obtained from the facets of the
            # simplicial complex that contain F by removing the vertices of F
            # from a facet that contains F.
            link_facets = []
            for facet in self.facets():
                if all(vertex in facet for vertex in F):
                    link_facets.append(tuple([v for v in facet if v not in F]))
            return SimplicialComplex(link_facets)

    def deletion(self, deleted_face):
        """Returns the simplicial complex obtained by deleting the given face
        from the simplicial complex. The original complex is unchanged.
        """

        # Form a new simplicial complex, deletion, whose faces are faces of
        # the original simplicial complex that do not contain the deleted face.
        deleted_complex = SimplicialComplex([])
        for i in range(-1, self.rank):
            for face in self.face_dictionary[i]:
                ##print(face, any(vertex not in face for vertex in deleted_face))
                # deleted face is not contained in face if there exists a
                # vertex in deleted_face that does not belong to face
                if any(vertex not in face for vertex in deleted_face):
                    deleted_complex.add_face(face)

        return deleted_complex

    def remove_face(self, deleted_face):
        """Removes the given face from the simplicial complex.
        """
        # Find all faces in the simplicial complex that contain the face that
        # is to be deleted. Remove them from the face dictionary.
        for i in range(-1, self.rank):
            # For each dimension, first locate the faces that need to be
            # deleted. They cannot be removed from the dictionary while
            # iterating over the set of elements of dimension i.
            faces_to_remove = set()
            for test_face in self.face_dictionary[i]:
                if all(vertex in test_face for vertex in deleted_face):
                    faces_to_remove.add(test_face)
            for face in faces_to_remove:
                self.face_dictionary[i].remove(face)

            # If deleting the given face results in the deletion of all faces
            # of dimension i, then remove that key from the face dictionary.
            if self.face_dictionary[i] == set():
                del self.face_dictionary[i]

        # If we deleted everything, add the empty face back in. We do not allow
        # for an empty complex.
        if len(self.face_dictionary.keys()) == 0:
            self.face_dictionary[-1] = set(tuple([]))

        # Update the dimension, rank, and vertex_set attributes of the
        # simplicial complex to account for the fact that we may have
        # changed the dimension or increased the number of vertices.
        self.dimension = max(self.face_dictionary.keys())
        self.rank = self.dimension + 1
        if self.dimension > -1:
            self.vertex_set = set([v[0] for v in self.face_dictionary[0]])
        else:
            self.vertex_set = set()
