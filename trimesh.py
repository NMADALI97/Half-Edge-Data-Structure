import numpy as np
from scipy.sparse import csr_matrix
from mpmath import *
import math
from scipy.sparse.linalg import eigsh

def get_heron_area(a, b, c):

    x = np.linalg.norm((b - a), 2)
    y = np.linalg.norm((c - a), 2)
    z = np.linalg.norm((c - b), 2)
    s = (x + y + z) * 0.5

    return (s * (s - x) * (s - y) * (s - z)) ** 0.5


    

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
class TriMesh( object ):
    def __init__( self ):
        self.edges = None
        self.vs = []
        self.faces = []


    def update_edge_list( self ):
        #from sets import Set, ImmutableSet
        Set, ImmutableSet = set, frozenset
        
        ## We need a set of set-pairs of vertices, because edges are bidirectional.
        edges = Set()
        for face in self.faces:
            edges.add( ImmutableSet( ( face[0], face[1] ) ) )
            edges.add( ImmutableSet( ( face[1], face[2] ) ) )
            edges.add( ImmutableSet( ( face[2], face[0] ) ) )
        
        self.edges = [ tuple( edge ) for edge in edges ]
    class HalfEdge( object ):
        def __init__( self ):
            self.to_vertex = -1
            self.face = -1
            self.edge = -1
            self.opposite_he = -1
            self.next_he = -1
    def update_halfedges( self ):
        
        self.halfedges = []
        self.vertex_halfedges = None
        self.face_halfedges = None
        self.edge_halfedges = None
        self.directed_edge2he_index = {}
        
        __directed_edge2face_index = {}
        for fi, face in enumerate( self.faces ):
            __directed_edge2face_index[ (face[0], face[1]) ] = fi
            __directed_edge2face_index[ (face[1], face[2]) ] = fi
            __directed_edge2face_index[ (face[2], face[0]) ] = fi
        
        def directed_edge2face_index( edge ):
            result = __directed_edge2face_index.get( edge, -1 )
            if -1 == result:
                assert edge[::-1] in __directed_edge2face_index
            
            return result
        
        self.vertex_halfedges = [None] * len( self.vs )
        self.face_halfedges = [None] * len( self.faces )
        self.edge_halfedges = [None] * len( self.edges )
        
        for ei, edge in enumerate( self.edges ):
            he0 = self.HalfEdge()
            
            he0.face = directed_edge2face_index( edge )
            he0.to_vertex = edge[1]
            he0.edge = ei
            
            he1 = self.HalfEdge()
           
            he1.face = directed_edge2face_index( edge[::-1] )
            he1.to_vertex = edge[0]
            he1.edge = ei
            
            
            he0index = len( self.halfedges )
            self.halfedges.append( he0 )
            he1index = len( self.halfedges )
            self.halfedges.append( he1 )
            
            
            he0.opposite_he = he1index
            he1.opposite_he = he0index
            
            
            assert edge not in self.directed_edge2he_index
            assert edge[::-1] not in self.directed_edge2he_index
            self.directed_edge2he_index[ edge ] = he0index
            self.directed_edge2he_index[ edge[::-1] ] = he1index
            
        
            if self.vertex_halfedges[ he0.to_vertex ] is None or -1 == he1.face:
                self.vertex_halfedges[ he0.to_vertex ] = he0.opposite_he
            if self.vertex_halfedges[ he1.to_vertex ] is None or -1 == he0.face:
                self.vertex_halfedges[ he1.to_vertex ] = he1.opposite_he
            
            
            if -1 != he0.face and self.face_halfedges[ he0.face ] is None:
                self.face_halfedges[ he0.face ] = he0index
            if -1 != he1.face and self.face_halfedges[ he1.face ] is None:
                self.face_halfedges[ he1.face ] = he1index
            
            
            assert self.edge_halfedges[ ei ] is None
            self.edge_halfedges[ ei ] = he0index
        
    
        boundary_heis = []
        for hei, he in enumerate( self.halfedges ):
            
            if -1 == he.face:
                boundary_heis.append( hei )
                continue
            
            face = self.faces[ he.face ]
            i = he.to_vertex
            j = face[ ( list(face).index( i ) + 1 ) % 3 ]
            
            he.next_he = self.directed_edge2he_index[ (i,j) ]
        

        vertex2outgoing_boundary_hei = {}
       
        Set = set
        for hei in boundary_heis:
            originating_vertex = self.halfedges[ self.halfedges[ hei ].opposite_he ].to_vertex
            vertex2outgoing_boundary_hei.setdefault(originating_vertex, Set()).add( hei )
            if len( vertex2outgoing_boundary_hei[ originating_vertex ] ) > 1:
                print('Butterfly vertex encountered')
        
        
        for hei in boundary_heis:
            he = self.halfedges[ hei ]
            for outgoing_hei in vertex2outgoing_boundary_hei[ he.to_vertex ]:
                he.next_he = outgoing_hei
                vertex2outgoing_boundary_hei[ he.to_vertex ].remove( outgoing_hei )
                break
        
        assert False not in [ 0 == len( out_heis ) for out_heis in vertex2outgoing_boundary_hei.values() ]
    def update_face_normals( self ):
    
        self.face_normals = np.zeros( ( len( self.faces ), 3 ) )
        self.face_areas = np.zeros( len( self.faces ) )
        
        vs = np.asarray( self.vs )
        fs = np.asarray( self.faces, dtype = int )
        
        
        self.face_normals = np.cross( vs[ fs[:,1] ] - vs[ fs[:,0] ], vs[ fs[:,2] ] - vs[ fs[:,1] ] )
        self.face_areas = np.sqrt((self.face_normals**2).sum(axis=1))
        self.face_normals /= self.face_areas[:,np.newaxis]
        self.face_areas *= 0.5

    def update_vertex_normals( self ):
        
        self.vertex_normals = np.zeros( ( len(self.vs), 3 ) )
        
        for vi in range( len( self.vs ) ):
            self.vertex_normals[vi] = 0.
            
            for fi in self.vertex_face_neighbors( vi ):
                
                self.vertex_normals[vi] += self.face_normals[ fi ] * self.face_areas[ fi ]
        
        
        self.vertex_normals *= 1./np.sqrt( ( self.vertex_normals**2 ).sum(1) ).reshape( (len(self.vs), 1) )
        
    def vertex_face_neighbors( self, vertex_index ):
        halfedges = self.halfedges
        result = []
        start_he = halfedges[ self.vertex_halfedges[ vertex_index ] ]
        he = start_he
        while True:
            if -1 != he.face: 
              result.append( he.face )
            
            he = halfedges[ halfedges[ he.opposite_he ].next_he ]
            if he is start_he: 
              break
        
        return result

    def vertex_vertex_neighbors( self, vertex_index ):
        
        halfedges = self.halfedges
        result = []
        start_he = halfedges[ self.vertex_halfedges[ vertex_index ] ]
        he = start_he
        while True:
            result.append( he.to_vertex )
            
            he = halfedges[ halfedges[ he.opposite_he ].next_he ]
            if he is start_he: 
              break
        
        return result
    
    
    def update_Laplacian_weights( self ):
 
        numv = self.vs.shape[0]
        numt = self.faces.shape[0]

        self.A = np.zeros((numv, numt))

        self.L = np.zeros((numv, numt, 3))

        for i in range(numv):

            req_t = self.faces[(self.faces[:, 0] == i) | (self.faces[:, 1] == i) | (self.faces[:, 2] == i)]

            for j in range(len(req_t)):

                tid = np.where(np.all(self.faces == req_t[j], axis=1))

                nbhr = [v for v in req_t[j] if v != i]

                vec1 = (self.vs[nbhr[0]] - self.vs[i]) /  np.linalg.norm(self.vs[nbhr[0]] - self.vs[i], 2)
                vec2 = (self.vs[nbhr[1]] - self.vs[i]) / np.linalg.norm(self.vs[nbhr[1]] - self.vs[i], 2)
                angle_at_x = np.arccos(np.dot(vec1, vec2))

                if angle_at_x > np.pi / 2:
                    self.A[i, tid] = get_heron_area(
                        self.vs[i], self.vs[nbhr[0]], self.vs[nbhr[1]]) / 2
                    continue

                vec1a = (self.vs[i] - self.vs[nbhr[0]]) / np.linalg.norm(self.vs[i] - self.vs[nbhr[0]], 2)
                vec2a = (self.vs[nbhr[1]] - self.vs[nbhr[0]]) / np.linalg.norm(self.vs[nbhr[1]] - self.vs[nbhr[0]], 2)

                inner_prod = np.dot(vec1a, vec2a)
                angle1 = np.arccos(inner_prod)

                if angle1 > np.pi / 2:
                    self.A[i, tid] = get_heron_area(self.vs[i], self.vs[nbhr[0]], self.vs[nbhr[1]]) / 4
                    continue

                vec1b = (self.vs[i] - self.vs[nbhr[1]]) / np.linalg.norm(self.vs[i] - self.vs[nbhr[1]], 2)
                vec2b = (self.vs[nbhr[0]] - self.vs[nbhr[1]]) / np.linalg.norm(self.vs[nbhr[0]] - self.vs[nbhr[1]], 2)

                inner_prod = np.dot(vec1b, vec2b)
                angle2 = np.arccos(inner_prod)

                if angle2 > np.pi / 2:
                    self.A[i, tid] = get_heron_area(
                        self.vs[i], self.vs[nbhr[0]], self.vs[nbhr[1]]) / 4
                    continue

                cot_1 = 1 / np.tan(angle1)
                cot_2 = 1 / np.tan(angle2)

                A_v_of_tid = 0.125 * ((cot_1 * np.linalg.norm(self.vs[i] - self.vs[nbhr[1]], 2)**2) + (cot_2 * np.linalg.norm(self.vs[i] - self.vs[nbhr[0]], 2)**2))

                self.L_at_v_t = ((1 / np.tan(angle1)) * (self.vs[i] - self.vs[nbhr[1]])) + ((1 / np.tan(angle2)) * (self.vs[i] - self.vs[nbhr[0]]))

                self.A[i, tid] = A_v_of_tid
                self.L[i, tid] = self.L_at_v_t

        self.A = np.sum(self.A, axis=1)
        # Set zeros in self.A to very small values
        self.A[self.A == 0] = 10 ** -40
        self.L = ((1 / (2 * self.A)) * np.sum(self.L, axis=1).T).T

        self. K_H = 0.5 * np.linalg.norm(self.L, 2, axis=1)

        
    
    
        
        
def TriMesh_From_OFF( path ):
        
    mesh= TriMesh()
    file = open( path )

    n_verts, n_faces, n_dontknow = tuple(
        [int(s) for s in file.readline().strip().split()])
    verts = [[float(s) for s in file.readline().strip().split()]
             for i_vert in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split()][1:]
             for i_face in range(n_faces)]
    
    mesh.vs=np.array(verts)

    mesh.faces=np.array(faces,dtype=int)
    mesh.update_edge_list( )
    mesh.update_halfedges()
    mesh.update_face_normals()
    mesh.update_vertex_normals()
    return mesh



