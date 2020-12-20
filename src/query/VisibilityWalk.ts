import { Vector } from "../types";
import { Facet } from "../geom/Geometry";

/**
 * Provides an initial estimate to start searching, based on the facets axis-oriented bounds.
 *
 * @ignore
 */
function findStartFacet(position: Vector, points: Vector[], facets: Facet[]): number
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

        if (found) return i;
    }

    return -1;
}

/**
 * Walks recursively through the neihbors of a set until the containing facet is found
 *
 * @ignore
 */
function walk(position: Vector, facet: Facet, points: Vector[]): Facet
{
    let foundRidge;

    for (let r of facet.ridges) {

    }

    if (!foundRidge)
        return facet;
    else if (!foundRidge.neighbor)
        return null;
    else
        return walk(position, foundRidge.neighbor.facet, points);
}


/**
 * Performs the visibility walk algorithm to find the Facet containing the given position. This should only be used
 * on Delaunay triangulations, as other triangulations are not guaranteed to resolve to a solution.
 * @param position The position to search for.
 * @param facets The facets to search
 * @param points The points indexed by the facets.
 * @param startIndex An optional index into the facets to start the search. If -1 is provided, an initial search
 * estimate may be made, but this is not guaranteed to be a performance improvement.
 */
export function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], startIndex: number = 0): Facet
{
    if (startIndex === -1) {
        startIndex = findStartFacet(position, points, facets);
        if (startIndex === -1)
            return null;
    }

    return walk(position, facets[startIndex], points);
}