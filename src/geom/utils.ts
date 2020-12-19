import { Facet, Ridge } from "./Geometry";
import { Vector } from "../types";
import { hyperplaneFromPoints, negate, signedDistToPlane } from "../math/VecMath";

/**
 * Some shape construction code used internally.
 *
 * @author derschmale <http://www.derschmale.com>
 */

/**
 * Assign the neighbor from a set. Generally only used while constructing a new facet.
 *
 * @param facet The facet owning the ridge.
 * @param ridge The ridge for which to find the neighbor.
 * @param facets The set of facets to search.
 *
 * @ignore
 */
export function findNeighbor(facet: Facet, ridge: Ridge, facets: Facet[]): void
{
    const src = ridge.verts;
    const numVerts = src.length;

    for (let f of facets) {
        for (let r of f.ridges) {
            // do not bother if we already found them
            if (r.neighbor) continue;

            let found = true;

            for (let i = 0; i < numVerts; ++i) {
                found = found && r.verts.indexOf(src[i]) >= 0;
            }

            // all vertices are shared
            if (found) {
                ridge.neighbor = r;
                r.neighbor = ridge;
                return;
            }
        }
    }
}

/**
 * Generates a plane for a facet. Used when constructing new facets.
 * @param facet The facet to construct a new plane for.
 * @param points The general set of points the indices refer to.
 * @param centroid An optional point that's guaranteed to be behind the plane. This can be used to orient the face
 * if the vertex order is inconsistent.
 *
 * @ignore
 */
export function generateFacetPlane(facet: Facet, points: Vector[], dim: number, centroid?: Vector): void
{
    const verts = facet.ridges.map(r => points[r.verts[0]]);

    let plane = facet.plane = hyperplaneFromPoints(verts);


    if (centroid && signedDistToPlane(centroid, plane, dim) > 0.0) {
        negate(plane);

        // flip ridges for consistency
        facet.ridges.reverse();

        for (let r of facet.ridges)
            r.verts.reverse();
    }
}