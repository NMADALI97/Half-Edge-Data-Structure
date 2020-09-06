import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
from trimesh import *



def main():



    
    trimesh=TriMesh_From_OFF( "torus.off")
    trimesh.update_Laplacian_weights()

    x, y, z = trimesh.vs[:, 0], trimesh.vs[:, 1], trimesh.vs[:, 2]

    mesh = mlab.triangular_mesh(x, y, z, trimesh.faces,
                                    representation='wireframe', opacity=0)

    mesh.mlab_source.dataset.point_data.scalars =trimesh.K_H
    mesh.mlab_source.dataset.point_data.scalars.name = 'Point data'

    mesh.mlab_source.update()
    mesh.parent.update()

    mesh2 = mlab.pipeline.set_active_attribute(mesh,
                                                   point_scalars='Point data')
    s2 = mlab.pipeline.surface(mesh2)
    s2.actor.mapper.interpolate_scalars_before_mapping = True
    mlab.colorbar(s2, title='Curvature\n', orientation='vertical')
    mlab.show()


    
main()