# Tympanum

A Typescript library for multidimensional geometry code. There's not much in it, but may grow slightly as I need more
functionality.

# Examples:

- [2D Quickhull](https://derschmale.github.io/tympanum/examples/convex_hull_2d/index.html)
- [3D Quickhull](https://derschmale.github.io/tympanum/examples/convex_hull_3d/index.html)

## Basic Types

Tympanum has the following building blocks to form shapes:

Any N-dimensional shape such as a simplex is a collection of Facets.
- `Facet`: This is a polygonal face of dimension N-1: a line, a triangle, or a tetrahedron in 2D, 3D or 4D respectively. 
  Each facet is bounded by a set of ridges.
- `Ridge`: This is an edge of dimension N-2: a point (vertex), a line (edge), or a triangle in 2D, 3D or 4D respectively.
  A ridge has N-1 vertices (ie: 1 vertex, 2 line end points, 3 triangle corners).
- `Vertex`: These are represented as an index into a list of points (fe: the list of points used to generate a convex 
  hull). This is so that we can easily map points to other data sets from which the points were extracted, or they can
  be used to construct 3D meshes for use in WebGL.
  
## Convex Hull

To generate a convex hull using the quickHull algorithm:

```
import { quickHull } from "@derschmale/tympanum";

let points = [];

for (let i = 0; i < 5000; ++i) {  
    points[i] = [Math.random(), Math.random(), Math.random()];
}

let hull = quickHull(points);

```

`hull` will contain an array of `Facet`.
