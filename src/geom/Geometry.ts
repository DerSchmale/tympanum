import { Vector } from "../types";

/**
 * Base geometry elements.
 *
 * @author derschmale <http://www.derschmale.com>
 */


/**
 * Ridge is a face of a facet, represented as the set of vertices forming the ridge. A line facet contains 2 ridges
 * each containing a single point. A triangle facet contains 3 ridges containing 2 vertices (segment end points), a
 * tetrahedral facet contains 4 face ridges (each containing the three triangle vertices).
 */
export class Ridge
{
    /**
     * The vertices of the ridge. These are always the points forming the ridge. Represented as indices into an external
     * point array.
     */
    verts: number[] = [];

    /**
     * The neighbour ridge belong to a different face but sharing the same vertices. In a closed shape, every ridge
     * should have a neighbour.
     */
    neighbor: Ridge;

    /**
     * The index of the opposite vertex. This is required by some query code.
     */
    opposite: number;

    /**
     * The facet to which this ridge belongs.
     */
    facet: Facet;

    /**
     * Creates a new ridge belonging to a facet
     */
    constructor(facet: Facet)
    {
        this.facet = facet;
    }
}

/**
 * In N dimensions, a facet forms an N-1 "polygon" which can be combined into an N-dimensional shape such as a
 * simplex, a convex hull, a triangulation, ...
 * A facet consists out of ridges and vertices. Care has to be taken during construction that the vertices and
 * ridges are in consistent order. The ridge's vertices must be a looping slice of (size-1) of the facet's vertices
 * starting at its corresponding index, ie:
 * - Every ridge at index I must start with the corresponding vertex at index I.
 * - If this ridge at index I stores a vertex at index N, the facet's vertex index must be (N + I) % numVerts.
 */
export class Facet
{
    /**
     * The set of ridges for the facet. Ridges are a dimension lower than the facet (ie: points for lines, edges for
     * triangles, faces for tetrahedrons).
     */
    ridges: Ridge[] = [];

    /**
     * The (hyper)plane containing the facet, represented as an N+1-dimensional vector (normal, offset)
     */
    plane: Vector;

    /**
     * The vertices of the facet. Represented as indices into an external point array. While they're also contained
     * in the ridges, it's useful for calculating barycentric coordinates for a facet.
     */
    verts: number[] = [];

    /**
     * Any sort of meta-data, generally used internally.
     */
    meta: any;
}