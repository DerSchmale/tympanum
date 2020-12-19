import { Vector } from "../types";
import { dim, hyperplaneFromPoints, signedDistToPlane } from "../math/VecMath";
import { Facet, Ridge } from "../geom/Geometry";
import { createSimplex } from "../geom/Simplex";
import { findNeighbor, generateFacetPlane } from "../geom/utils";
import { removeElementOutOfOrder, removeIndexOutOfOrder, shuffle } from "@derschmale/array-utils";

const eps = 0.0001;

/**
 * Some meta-data while constructing the facets
 */
class FacetInfo
{
    outsideSet: number[] = [];
    outsideDist: number[] = [];
    currentPoint: Vector;   // used to keep track of visibility tests
}

/**
 * Gets the point in an outside set that's the furthest from the facet plane.
 *
 * @ignore
 */
function getFurthestPoint(facet: Facet): number
{
    const { outsideSet, outsideDist } = facet.meta;
    const len = outsideSet.length;

    if (len === 0) return -1;

    let p = outsideSet[0];
    let dist = outsideDist[0];

    for (let i = 1; i < len; ++i) {
        if (outsideDist[i] > dist) {
            dist = outsideDist[i];
            p = outsideSet[i];
        }
    }

    return p;
}

/**
 * Assigns all points to the outside sets of a face.
 *
 * @ignore
 */
function generateOutsideSets(indices: number[], points: Vector[], facets: Facet[], dim: number)
{
    let outsideSets = facets.map(_ => []);

    const len = indices.length;
    for (let i = 0; i < len; ++i) {
        let index = indices[i];
        const p = points[index];

        for (let f of facets) {
            let dist = signedDistToPlane(p, f.plane, dim);

            if (dist > eps) {
                const meta = f.meta;
                meta.outsideSet.push(index);
                meta.outsideDist.push(dist);
                break;
            }
        }
    }

    // any point now not in an outside set is inside the hull
    return outsideSets;
}

/**
 * Finds all faces visible to a point and their boundary ridges.
 *
 * @ignore
 */
function getVisibleSet(p: Vector, facet: Facet, visible: Facet[], horizon: Ridge[], dim: number)
{
    visible.push(facet);
    facet.meta.currentPoint = p;

    for (let r of facet.ridges) {
        const neighbor = r.neighbor.facet;

        // already checked
        if (neighbor.meta.currentPoint === p) continue;

        if (signedDistToPlane(p, neighbor.plane, dim) > eps)
            getVisibleSet(p, neighbor, visible, horizon, dim);
        else
            horizon.push(r);
    }
}

/**
 * Builds a new face from a ridge and a point.
 *
 * @ignore
 */
function attachNewFacet(ridge: Ridge, p: number, points: Vector[], facets: Facet[], centroid: Vector, dim: number): Facet
{
    // in 2D, we simply need to create 1 new facet (line) from old ridge to p
    const newFacet = new Facet();
    newFacet.meta = new FacetInfo();

    // collect all verts for this facet, which is the horizon ridge + this point
    const verts = ridge.verts.slice();
    verts.push(p);

    // the horizon ridge is part of the new facet, and gets to keep its neighbor
    newFacet.ridges.push(ridge);
    ridge.facet = newFacet;

    // dim + 1 ridges (3 edges to a triangle in 2D, 4 faces to a tetrahedron in 3D)
    // start with 1, since 0 would be the same as the already existing ridge above
    for (let r = 1; r < dim; ++r) {
        let ridge = new Ridge(newFacet);

        for (let v = 0; v < dim - 1; ++v)
            ridge.verts[v] = verts[(r + v) % dim];

        // so we only need to search for neighbours in the newly generated facets, only the horizons are attached to
        // the old ones
        findNeighbor(newFacet, ridge, facets);

        newFacet.ridges.push(ridge);
    }

    generateFacetPlane(newFacet, points, dim, centroid);

    return newFacet;
}

/**
 * Builds a set of new facets for a point and its horizon ridges.
 *
 * @ignore
 */
function connectHorizonRidges(points: Vector[], index: number, H: Ridge[], centroid: Vector, dim: number): Facet[]
{
    const newFacets = [];

    // link horizon ridges with new point
    for (let ridge of H) {
        const newFacet = attachNewFacet(ridge, index, points, newFacets, centroid, dim);
        newFacets.push(newFacet);
    }

    return newFacets;
}

/**
 * Creates the centroid for a collection of d points.
 *
 * @ignore
 */
function createCentroid(points: Vector[], d: number, indices: number[]): Vector
{
    // a point that will be internal from the very first simplex. Used to correctly orient new planes
    const centroid = points[0].slice();

    for (let j = 0; j < d; ++j) {
        for (let i = 1; i <= d; ++i) {
            let index = indices[i];
            centroid[j] += points[index][j];
        }
        centroid[j] /= d + 1;
    }

    return centroid;
}

/**
 * Returns the index with the "largest" point. The largest point is the one with the highest x coefficient, or y, z,
 * etc. if equal
 *
 * @ignore
 */
function maximize(i: number, maxIndex: number, points: Vector[], d: number): number
{
    const p = points[i];
    const max = points[maxIndex];

    for (let j = 0; j < d; ++j) {
        if (p[j] < max[j]) return maxIndex;
        if (p[j] > max[j]) return i;
    }

    // all the same: this only happens with duplicates, which shouldn't be in the set
    return maxIndex;
}

/**
 * @ignore
 */
function minimize(i: number, minIndex: number, points: Vector[], d: number): number
{
    const p = points[i];
    const min = points[minIndex];

    for (let j = 0; j < d; ++j) {
        if (p[j] > min[j]) return minIndex;
        if (p[j] < min[j]) return i;
    }

    // all the same: this only happens with duplicates, which shouldn't be in the set
    return minIndex;
}

/**
 * Tries to find the biggest shape to start with
 * @ignore
 */
function getOptimalStart(points: Vector[], d: number): number[]
{
    const numPoints = points.length;
    let minIndex = 0;
    let maxIndex = 0;

    // the initial axis
    for (let i = 1; i < numPoints; ++i) {
        maxIndex = maximize(i, maxIndex, points, d);
        minIndex = minimize(i, minIndex, points, d);
    }

    const indices = [ minIndex, maxIndex ];
    const planePts = [ points[minIndex], points[maxIndex] ];

    // already have 2 points, need d + 1 in total
    // in increasing dimensions, find the furthest from the current hyperplane
    for (let i = 2; i < d + 1; ++i) {
        const plane = hyperplaneFromPoints(planePts);
        let maxDist = -Infinity;
        let p = -1;

        for (let j = 0; j < numPoints; ++j) {
            const dist = Math.abs(signedDistToPlane(points[j], plane, i));
            if (dist > maxDist) {
                maxDist = dist;
                p = j;
            }
        }

        indices.push(p);
        planePts.push(points[p]);
    }

    return indices;
}

/**
 * QuickHull implements the algorithm of the same name, based on the original paper by Barber, Dobkin and Huhdanpaa.
 * We're not interested in 0- or 1-dimensional cases (the latter can simply be sorted points). QuickHull returns a
 * set of indices into the original point list so we can map it to a different original ata set (fe: points may be a
 * mapping for position vectors on some scene graph object).
 *
 * @author derschmale <http://www.derschmale.com>
 */
export function quickHull(points: Vector[]): Facet[]
{
    const numPoints = points.length;
    if (numPoints === 0) return;
    const d = dim(points[0]);

    if (numPoints <= d) {
        console.log(`A convex hull in ${d} dimensions requires at least ${d + 1} points.`);
    }

    // initial unprocessed point indices:
    let indices = [];
    for (let i = 0; i < numPoints; ++i)
        indices.push(i);

    const simplexIndices = getOptimalStart(points, d);

    const centroid = createCentroid(points, d, simplexIndices);
    const facets = createSimplex(points, d, simplexIndices);

    for (let f of facets)
        f.meta = new FacetInfo();

    // sorting them in descending order makes them easy to delete optimally
    simplexIndices.sort((a, b) => b - a);
    for (let i = 0; i < d + 1; ++i) {
        removeIndexOutOfOrder(indices, simplexIndices[i]);
    }

    shuffle(indices);

    generateOutsideSets(indices, points, facets, d);

    // do not cache facets.length, as it will keep changing
    let done = false;

    // TODO: this extra loop should not be required
    while (!done) {
        done = true;
        for (let i = 0; i < facets.length; ++i) {
            const facet = facets[i];

            const p = getFurthestPoint(facet);
            if (p !== -1) {
                removeElementOutOfOrder(facet.meta.outsideSet, p);

                const V = [];
                const H = [];

                getVisibleSet(points[p], facet, V, H, d);

                const newFacets = connectHorizonRidges(points, p, H, centroid, d);

                for (let v of V) {
                    if (removeElementOutOfOrder(facets, v) <= i) --i;
                    generateOutsideSets(v.meta.outsideSet, points, newFacets, d);
                    if (v.meta.outsideSet.length > 0)
                        done = false;
                }

                facets.push(...newFacets);
            }
        }
    }


    for (let f of facets)
        f.meta = null;

    return facets;
}