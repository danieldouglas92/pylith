# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.testing.FullTestApp import MeshEntity

class TriGmsh(object):
    """Mesh information for tri mesh generated by Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=184, ncorners=3, nvertices=81),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "elastic_xneg": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "elastic_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh generated by Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=99, ncorners=4, nvertices=120),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        "elastic_xneg": MeshEntity(ncells=54, ncorners=4, nvertices=70),
        "elastic_xpos": MeshEntity(ncells=45, ncorners=4, nvertices=60),

        "bc_xneg": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_xpos": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_yneg": MeshEntity(ncells=11, ncorners=2, nvertices=12),
        "bc_ypos": MeshEntity(ncells=11, ncorners=2, nvertices=12),
    }


class TriCubit(object):
    """Mesh information for tri mesh generated by Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=124, ncorners=3, nvertices=79),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "elastic_xneg": MeshEntity(ncells=60, ncorners=3, nvertices=43),
        "elastic_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class QuadCubit(object):
    """Mesh information for quad mesh generated by Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=4, nvertices=81),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        "elastic_xneg": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "elastic_xpos": MeshEntity(ncells=32, ncorners=4, nvertices=45),

        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }



# End of file
