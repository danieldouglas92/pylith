#!/usr/bin/env python3
"""Generate a tri or quad mesh of a strike-slip fault using Gmsh, making use of the
built-in geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.

Run `generate_gmsh.py --write` to generate the mesh.
"""

# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.

    Domain is 100km by 150km.
    -50.0 km <= x <= 50.0 km
    -75.0 km <= y <= 75.0 km

    The fault surface runs along the y-axis through the entire domain.

    p4-----p6-----p3
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    p1-----p5-----p2
    """
    DOMAIN_X = 800.0e+3
    DOMAIN_Y = 400.0e+3
    LAYER_THICKNESS = 50.0e+3


    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "my_mesh.msh"

    def create_geometry(self):
        """Create geometry.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set local variables for domain size and corner of the domain.
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x1 = -0.5 * lx
        y1 = -1.0 * ly

        # Create points.
        self.p1 = gmsh.model.geo.add_point(x1, y1, 0.0)
        self.p2 = gmsh.model.geo.add_point(x1+lx, y1, 0.0)
        self.p3 = gmsh.model.geo.add_point(x1+lx, y1+ly, 0.0)
        self.p4 = gmsh.model.geo.add_point(x1, y1+ly, 0.0)
        self.p5 = gmsh.model.geo.add_point(x1, -self.LAYER_THICKNESS, 0.0)
        self.p6 = gmsh.model.geo.add_point(x1+lx, -self.LAYER_THICKNESS, 0.0)

        # Create Trapezoid points.
        self.p7 = gmsh.model.geo.add_point(0.5 * x1, y1 + ly, 0.0)
        self.p8 = gmsh.model.geo.add_point(0.25 * x1, y1 + 0.75 * ly, 0.0)
        self.p9 = gmsh.model.geo.add_point(-0.25 * x1, y1 + 0.75 * ly, 0.0)
        self.p10 = gmsh.model.geo.add_point(-0.5 * x1, y1 + ly, 0.0)
        self.p11 = gmsh.model.geo.add_point(0.375 * x1, -self.LAYER_THICKNESS, 0.0)
        self.p12 = gmsh.model.geo.add_point(-0.375 * x1, -self.LAYER_THICKNESS, 0.0)


        # p5 = gmsh.model.geo.add_point(x1+0.5*lx, y1, 0.0)
        # p6 = gmsh.model.geo.add_point(x1+0.5*lx, y1+ly, 0.0)

        # Create curves. We store the curve tag as a data member
        # so that we can refer to them later.
        self.c_bot = gmsh.model.geo.add_line(self.p1, self.p2)

        self.c_right_lower = gmsh.model.geo.add_line(self.p2, self.p6)
        self.c_right_layer = gmsh.model.geo.add_line(self.p6, self.p3)

        self.c_top_left = gmsh.model.geo.add_line(self.p7, self.p4)
        self.c_top_right = gmsh.model.geo.add_line(self.p3, self.p10)

        self.c_left_layer = gmsh.model.geo.add_line(self.p4, self.p5)
        self.c_left_lower = gmsh.model.geo.add_line(self.p5, self.p1)

        self.left_layer_bot = gmsh.model.geo.add_line(self.p11, self.p5)
        self.right_layer_bot = gmsh.model.geo.add_line(self.p12, self.p6)

        self.trap_top = gmsh.model.geo.add_line(self.p10, self.p7)

        self.trap_left_layer = gmsh.model.geo.add_line(self.p7, self.p11)
        self.trap_left_lower = gmsh.model.geo.add_line(self.p11, self.p8)

        self.trap_bot = gmsh.model.geo.add_line(self.p8, self.p9)

        self.trap_right_layer = gmsh.model.geo.add_line(self.p12, self. p10)
        self.trap_right_lower = gmsh.model.geo.add_line(self.p9, self.p12)

        # Create curve loops and surfaces from the curves.
        # We traverse the curves in a counter clock-wise direction.
        # If the curve is in the opporite direction, we use the negative tag.
        c_trap = gmsh.model.geo.add_curve_loop([self.trap_top, self.trap_left_layer, self.trap_left_lower, self.trap_bot, self.trap_right_lower, self.trap_right_layer])
        self.s_trap = gmsh.model.geo.add_plane_surface([c_trap])

        c_box = gmsh.model.geo.add_curve_loop([-self.trap_bot, -self.trap_left_lower, self.left_layer_bot, self.c_left_lower, self.c_bot, self.c_right_lower, -self.right_layer_bot, -self.trap_right_lower])
        self.s_box = gmsh.model.geo.add_plane_surface([c_box])
        
        c_layer_left = gmsh.model.geo.add_curve_loop([self.c_top_left, self.c_left_layer, -self.left_layer_bot, -self.trap_left_layer])
        self.s_left_layer = gmsh.model.geo.add_plane_surface([c_layer_left])

        c_layer_right = gmsh.model.geo.add_curve_loop([self.c_top_right, -self.trap_right_layer, self.right_layer_bot, self.c_right_layer])
        self.s_right_layer = gmsh.model.geo.add_plane_surface([c_layer_right])

        # curves_left = gmsh.model.geo.split_curve(self.c_left, [p5])
        # curves_right = gmsh.model.geo.split_curve(self.c_right, [p6])

        # self.c_topo_left_1 = curves_left[0]
        # self.c_topo_left_2 = curves_left[1]

        # self.c_topo_right_1 = curves_right[0]
        # self.c_topo_right_2 = curves_right[1]

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create two materials, one for each side of the fault.
        # We use the `MaterialGroup` data class defined in `gmsh_utils.`
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_box]),
            MaterialGroup(tag=2, entities=[self.s_trap]),
            MaterialGroup(tag=3, entities=[self.s_left_layer]),
            MaterialGroup(tag=4, entities=[self.s_right_layer])
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        # We use the `VertexGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.

        vertex_groups = (
            VertexGroup(name="boundary_ypos", tag=10, dim=1, entities=[self.c_top_left, self.trap_top, self.c_top_right]),
            VertexGroup(name="boundary_xneg", tag=11, dim=1, entities=[self.c_left_layer, self.c_left_lower]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_bot]),
            VertexGroup(name="boundary_xpos", tag=13, dim=1, entities=[self.c_right_layer, self.c_right_lower]),
            VertexGroup(name="fault", tag=20, dim=1, entities=[self.trap_left_layer, self.trap_left_lower]),
            VertexGroup(name="fault_edge", tag=21, dim=0, entities=[self.p8]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set discretization size with geometric progression from distance to the fault.
        
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.trap_left_layer, self.trap_left_lower])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=4.0e+3, bias=1.05)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()


# End of file
