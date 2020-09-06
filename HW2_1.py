import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

mappable = cm.ScalarMappable(cmap=cm.coolwarm, norm=matplotlib.colors.Normalize(vmin=0, vmax=1))

from trimesh import *
import glob
from PIL import Image
import pygame
from pygame.locals import *
import os
from OpenGL.GL import *
from OpenGL.GLU import *

mesh=TriMesh_From_OFF( "torus.off")

mesh.update_Laplacian_weights()
Colors = mappable.to_rgba(mesh.K_H)
print(Colors.shape)

def Curvature( ):
    glPolygonOffset( 1.0, 1.0 )
    glEnable( GL_POLYGON_OFFSET_FILL )
    
    glBegin( GL_TRIANGLES )
    for face, normal in zip( mesh.faces, mesh.face_normals ):
        glNormal3f( *normal )
        
        for vertex_index in face:
            
            glColor3f(*(Colors[vertex_index][1:]))
            glVertex3f( *mesh.vs[ vertex_index ] )
    glEnd()
    glDisable( GL_POLYGON_OFFSET_FILL )
def Mesh( ):
    glPolygonOffset( 1.0, 1.0 )
    glEnable( GL_POLYGON_OFFSET_FILL )
    
    glBegin( GL_TRIANGLES )
    for face, normal in zip( mesh.faces, mesh.face_normals ):
        glNormal3f( *normal )
        for vertex_index in face:
            glColor3f(0.50,0.50,0.50)
            glVertex3f( *mesh.vs[ vertex_index ] )
    glEnd()
    glDisable( GL_POLYGON_OFFSET_FILL )
    
    glBegin(GL_LINES)   
    for face, normal in zip( mesh.faces, mesh.face_normals ):     
        center=[0,0,0]
        for vertex_index in face:
            center+=mesh.vs[ vertex_index ]
        center/=3
        glColor3f(1.0, 0.0, 0.0)
        glVertex3f(*center)
        glVertex3f(*(center+normal))
    glEnd()
    #Surface Normal
    glBegin(GL_LINES)
    for edge in mesh.edges:
        for vertex in edge:
            glColor3f(1.0, 1.0, 1.0)
            glVertex3fv(mesh.vs[vertex])
    glEnd()
    #Vertex Normal
    glBegin(GL_LINES)
    for v, normal in zip( mesh.vs, mesh.vertex_normals ):     
        
        glColor3f(0.0, 1.0, 0.0)
        glVertex3f(*v)
        glVertex3f(*(v+normal))
    glEnd()
    
    
    

def main():
    cnt=0
    pygame.init()
    display = (800,600)
    screen = pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

    gluPerspective(45, (display[0]/display[1]), 0.1, 50.0)

    glTranslatef(0.0,0.0, -10)

    while True:
        for event in pygame.event.get():
            
            if event.type == pygame.QUIT:
               
               
                pygame.quit()
                quit()

        glRotatef(1, 3, 1, 1)
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        Curvature( )
        pygame.display.flip()
        pygame.time.wait(1)
       


main()