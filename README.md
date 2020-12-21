# Tympanum

A Typescript library for multidimensional geometry code. There's not much in it, but may grow slightly as I need more
functionality.

# Documentation

- [Reference documentation](https://derschmale.github.io/tympanum/docs/index.html)
- [The original QuickHull algorithm](https://www.researchgate.net/publication/2641780_The_QuickHull_Algorithm_for_Convex_Hulls)

# Examples:

- [2D Quickhull](https://derschmale.github.io/tympanum/examples/convex_hull_2d/index.html)
- [3D Quickhull](https://derschmale.github.io/tympanum/examples/convex_hull_3d/index.html)
- [2D Delaunay Triangulation](https://derschmale.github.io/tympanum/examples/delaunay_2d/index.html)
- [3D Delaunay Tetrahedralisation](https://derschmale.github.io/tympanum/examples/delaunay_3d/index.html)
- [Delaunay facet search using visibility walking](https://derschmale.github.io/tympanum/examples/walk_2d/index.html)
- [Point reconstruction with barycentric coordinates](https://derschmale.github.io/tympanum/examples/barycentric/index.html)

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

const points = [];

for (let i = 0; i < 5000; ++i) {  
    points[i] = [Math.random(), Math.random(), Math.random()];
}

const hull = quickHull(points);

```

`hull` will contain an array of `Facet`.

## Delaunay Triangulation

To generate the delaunay triangulation:

```
import { delaunay } from "@derschmale/tympanum";

const points = [];

for (let i = 0; i < 500; ++i) {  
    points[i] = [Math.random(), Math.random(), Math.random()];
}

const triangulation = delaunay(points);

```

`triangulation` will contain an array of `Facet`, but of a higher dimension than the convex hull would.

Delaunay triangulations allow searching for facets containing a point efficiently using the vibility walk algorithm:

```
import { visibilityWalk } from "@derschmale/tympanum";

const pos = [ 0.5, 0.2, 0.7 ];
const facet = visibilityWalk(pos, triangulation, points);

```

When a facet has been found, we can calculate the point's barycentric coordinates. The barycentric coordinates can be 
used to interpolate values associated to each respective point.

```
import { barycentricCoords } from "@derschmale/tympanum";

// for example: every point has an RGB color assigned to it:
let colors = [];

// any color at index N is associated with the point at points[N]
for (let i = 0; i < 5000; ++i) {  
    colors[i] = { 
      r: Math.random() * 0xff, 
      g: Math.random() * 0xff, 
      b: Math.random() * 0xff
    };
}

if (facet) {
  const bary = barycentricCoords(pos, facet, points);
  const color = { r: 0, g: 0, b: 0 };
  
  for (let i = 0; i < bary.length; ++i) {
    // get the index of the point
    let index = facet.verts[i];
  
    // get the color at that index
    let c = colors[index];
    
    // add the weighted colors together
    color.r += bary[i] * c.r; 
    color.g += bary[i] * c.g; 
    color.b += bary[i] * c.b; 
  }
}

```
