import { Facet, Ridge } from "./Geometry";
import { Vector } from "../types";
import { hyperplaneFromPoints, negate, signedDistToPlane } from "../math/VecMath";

/**
 * Some shape construction and query code used internally.
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
    // this simply checks if all vertices are shared for all provided facets. Not very efficient, there's probably
    // better ways to do a neighbour search. However, this is generally only applied to relatively small sets (newly
    // constructed faces).
    const src = ridge.verts;
    const numVerts = src.length;

    for (let f of facets) {
        for (let r of f.ridges) {
            // do not bother if we already found them
            if (r.neighbor) continue;

            let found = true;

            for (let i = 0; i < numVerts; ++i)
                found = found && r.verts.indexOf(src[i]) >= 0;

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

        facet.verts.reverse();

        // flip ridges for consistency
        facet.ridges.reverse();

        for (let r of facet.ridges)
            r.verts.reverse();
    }
}

/**
 * Generates the ridges for a facet based on the vertices contained in it. If ridges are already present, they're
 * considered already built and only the missing ones will be added (used only when connecting ridges). In this
 * case, the order of vertices MUST match that of the present ridges.
 *
 * @ignore
 */
export function buildRidges(facet: Facet, facets: Facet[], dim: number)
{
    const verts = facet.verts;
    const ridges = facet.ridges;

    for (let r = ridges.length; r < dim; ++r) {
        let ridge = new Ridge(facet);

        for (let v = 0; v < dim - 1; ++v) {
            ridge.verts[v] = verts[(r + v) % dim];
        }

        ridge.opposite = verts[(r + dim) % dim];
        findNeighbor(facet, ridge, facets);
        facet.ridges.push(ridge);
    }
}

/**
 * Combines a ridge and a point into a new facet.
 * @param ridge The ridge to extend.
 * @param p The index of the new point to build the missing ridges from.
 * @param points The array containing the point values.
 * @param facets The other facets in the shape, used to find neighbours. Usually, when constructing closed shapes,
 * these are only the newly constructed faces.
 * @param insidePoint A point guaranteed to be on the negative side of the facet plane.
 * @param dim The dimension of the facet.
 *
 * @ignore
 */
export function extendRidge(ridge: Ridge, p: number, points: Vector[], facets: Facet[], insidePoint: Vector, dim: number): Facet
{
    // in 2D, we simply need to create 1 new facet (line) from old ridge to p
    const facet = new Facet();

    // collect all verts for this facet, which is the horizon ridge + this point
    facet.verts = ridge.verts.concat([p]);

    // the horizon ridge is part of the new facet, and gets to keep its neighbor
    facet.ridges.push(ridge);
    ridge.facet = facet;

    // the new opposite point
    ridge.opposite = p;

    buildRidges(facet, facets, dim);
    generateFacetPlane(facet, points, dim, insidePoint);
    return facet;
}

/**
 * Calculates the centroid ("average") for a collection of points.
 *
 * @param points The array containing all point coordinates.
 * @param indices The indices of the points for which to calculate the averages. If not provided, the first N+1
 * (simplex) points are used from the points array.
 * @ignore
 */
export function createCentroid(points: Vector[], indices?: number[]): Vector
{
    // a point that will be internal from the very first simplex. Used to correctly orient new planes
    const centroid = points[0].slice();
    const dim = centroid.length;
    const numPoints = indices? indices.length : dim + 1;

    for (let j = 0; j < dim; ++j) {
        for (let i = 1; i < numPoints; ++i) {
            let index = indices[i];
            centroid[j] += points[index][j];
        }
        centroid[j] /= numPoints;
    }

    return centroid;
}