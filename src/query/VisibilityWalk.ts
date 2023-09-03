import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
import { dim, intersectRayPlane } from "../math/VecMath";
import { createCentroid } from "../geom/utils";
import { EPSILON } from "../constants";

/**
 * Provides an initial estimate to start searching, based on the facets axis-oriented bounds.
 *
 * @ignore
 */
function findStartFacet(position: Vector, points: Vector[], facets: Facet[]): Facet
{
    const numFacets = facets.length;
    for (let i = 0; i < numFacets; ++i) {
        const f = facets[i];
        const verts = f.verts;
        const numVerts = verts.length;
        const dim = points[0].length;
        let found = true;

        for (let d = 0; d < dim; ++d) {
            let min = points[verts[0]][d];
            let max = min;

            // find min coordinate
            for (let v = 1; v < numVerts; ++v) {
                const p = points[verts[v]][d];
                if (p < min) min = p;
                else if (p > max) max = p;
            }

            const p = position[d];

            if (p < min || p >= max) {
                found = false;
                break;
            }
        }

        if (found) return f;
    }

    return null;
}

/**
 * Walks recursively through the neihbors of a set until the containing facet is found.
 *
 * @ignore
 */
function walk(position: Vector, facet: Facet, points: Vector[], centroid: Vector, dir: Vector, dim: number): Facet
{
    createCentroid(points, facet.verts, centroid);

    // so now we need to test the ray centroid -> position against the ridges and see if any intersect
    for (let i = 0; i < dim; ++i) {
        dir[i] = position[i] - centroid[i];
    }

    // using the centroid makes things easier, as the ray starting from the centroid hits the triangle face for
    // which the intersection distance is closest, so test for minT rather than doing barycentric tests.
    let hit = null;
    let minT = 1.0 - EPSILON;

    for (let r of facet.ridges) {
        const t = intersectRayPlane(centroid, dir, r.getPlane(points, centroid), dim, true);

        // intersection did not occur on the segment, or it's not the furthest
        if (t > EPSILON && t <= minT) {
            minT = t;
            hit = r;
        }
    }

    // if no intersection is found, we must be in the facet
    if (hit) {
        // there is an intersection, but we may have left the shape if there's no neighbour
        return hit.neighbor ?
            walk(position, hit.neighbor.facet, points, centroid, dir, dim) :
            null;
    }

    return facet;
}


/**
 * Performs the visibility walk algorithm to find the Facet containing the given position. This should only be used
 * on Delaunay triangulations, as other triangulations are not guaranteed to resolve to a solution.
 * @param position The position to search for.
 * @param facets The facets to search
 * @param points The points indexed by the facets.
 * @param startFacet An optional facet to start the search. If none is provided, an initial search
 * estimate is made, but this is not guaranteed to be a performance improvement.
 */
export function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], startFacet?: Facet): Facet;


/**
 * Performs the visibility walk algorithm to find the Facet containing the given position. This should only be used
 * on Delaunay triangulations, as other triangulations are not guaranteed to resolve to a solution.
 * @param position The position to search for.
 * @param facets The facets to search
 * @param points The points indexed by the facets.
 * @param estimate If true, searches the facets to find an initial match.
 */
export function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], estimate?: boolean): Facet;

export function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], startFacetOrEstimate?: Facet | boolean): Facet
{
    let startFacet;
    if (startFacetOrEstimate === true) {
        startFacet = findStartFacet(position, points, facets);

        if (!startFacet)
            return null;
    }
    else startFacet = startFacetOrEstimate || facets[0];

    // this is just a single reusable object in order not to have to recreate it
    const d = dim(points[0]);
    const centroid = new Float32Array(d);
    const dir = new Float32Array(d);
    return walk(position, startFacet, points, centroid, dir, d);
}
