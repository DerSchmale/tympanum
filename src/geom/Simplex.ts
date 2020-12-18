import { Vector } from "../types";
import { Facet, Ridge } from "./Geometry";
import { findNeighbor, generateFacetPlane } from "./utils";

/**
 * Creates an N-simplex from N+1 points.
 *
 * @param points An array of points. Only the first N+1 points are used
 * @param dim The dimension of the simplex
 */
export function createSimplex(points: Vector[], dim: number): Facet[]
{
    // TODO: We can find a better initial data set, similar to QHull:
    //  find the minX, maxX points, and iteratively extend with furthest point
    //  (starts with signed dist to line, then to plane in 3D, then to hyperplane in 4D)
    //  This is the same sort of logic of the base Quickhull algorithm, so maybe it's not that much of an
    //  improvement to do it in the first step?


    const facets: Facet[] = [];
    const numVerts = dim + 1;
    const verts = [];

    for (let i = 0; i <= dim; ++i) {
        const f = new Facet();

        // collect all verts for this facet
        // the facet is made up of dim + 1 points, so cycle through these
        for (let v = 0; v < dim; ++v) {
            verts[v] = (i + v) % numVerts;
        }

        for (let r = 0; r < dim; ++r) {

            let ridge = f.ridges[r] = new Ridge(f);

            for (let v = 0; v < dim - 1; ++v) {
                ridge.verts[v] = verts[(r + v) % dim];
            }

            // find neighbours from already generated sets
            // there's probably an analytical way to do this
            findNeighbor(f, ridge, facets);
        }

        // an opposing face to ensure correct direction
        const opposing = points[(i + dim) % (dim + 1)];
        generateFacetPlane(f, points, opposing);

        facets.push(f);
    }

    return facets;
}