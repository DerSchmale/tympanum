import { Vector } from "../types";
import { Facet } from "./Geometry";
import { buildRidges, generateFacetPlane } from "./utils";

/**
 * Creates an N-simplex from N+1 points. The dimension of the points are used to
 * define the dimension of the simplex.
 *
 * @param points An array of points. Only the first N+1 points are used.
 * @param indices An optional array of indices into points to define which points in the set are used.
 *
 * @author derschmale <http://www.derschmale.com>
 */
export function createSimplex(points: Vector[], indices?: number[]): Facet[]
{
    const dim = points[0].length;
    const numVerts = dim + 1;
    const facets: Facet[] = [];

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