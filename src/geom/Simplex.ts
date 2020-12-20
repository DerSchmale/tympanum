import { Vector } from "../types";
import { Facet } from "./Geometry";
import { buildRidges, generateFacetPlane } from "./utils";

/**
 * Creates an N-simplex from N+1 points.
 *
 * @param points An array of points. Only the first N+1 points are used.
 * @param dim The dimension of the simplex.
 * @param indices An optional array of indices into points to remap which points are used
 *
 * @author derschmale <http://www.derschmale.com>
 */
export function createSimplex(points: Vector[], dim: number, indices?: number[]): Facet[]
{
    // TODO: We can find a better initial data set, similar to QHull:
    //  find the minX, maxX points, and iteratively extend with furthest point
    //  (starts with signed dist to line, then to plane in 3D, then to hyperplane in 4D)
    //  This is the same sort of logic of the base Quickhull algorithm, so maybe it's not that much of an
    //  improvement to do it in the first step?


    const facets: Facet[] = [];
    const numVerts = dim + 1;

    for (let i = 0; i <= dim; ++i) {
        const f = new Facet();

        // collect all verts for this facet
        // the facet is made up of dim + 1 points, so cycle through these
        for (let v = 0; v < dim; ++v) {
            const index = (i + v) % numVerts;
            f.verts[v] = indices ? indices[index] : index;
        }

        buildRidges(f, facets, dim);

        // an opposing face to ensure correct direction
        const index = (i + dim) % (dim + 1);
        const opposing = points[indices ? indices[index] : index];
        generateFacetPlane(f, points, dim, opposing);

        facets.push(f);
    }

    return facets;
}