import { Vector } from "../types";
import { dim, signedDistToPlane } from "../math/VecMath";
import { removeElementOutOfOrder } from "../array/utils";
import { Facet, Ridge } from "../geom/Geometry";
import { createSimplex } from "../geom/Simplex";
import { findNeighbor, generateFacetPlane } from "../geom/utils";

/**
 * Some meta-data while constructing the facets
 */
class FacetInfo
{
    outsideSet: number[] = [];
    outsideDist: number[] = [];
    currentPoint: Vector;   // used to keep track of visibility tests
}

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

// [face][setIndex][0/1] : 0 = index, 1 = signed distance to facet plane
// assign points to the outside set of a collection of faces
function generateOutsideSets(indices: number[], points: Vector[], facets: Facet[])
{
    let outsideSets = facets.map(_ => []);

    const len = indices.length;
    for (let i = 0; i < len; ++i) {
        let index = indices[i];
        const p = points[index];

        for (let f of facets) {
            let dist = signedDistToPlane(p, f.plane);

            if (dist > 0) {
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

function getVisibleSet(p: Vector, facet: Facet, visible: Facet[], horizon: Ridge[])
{
    visible.push(facet);
    facet.meta.currentPoint = p;

    for (let r of facet.ridges) {
        const neighbor = r.neighbor.facet;

        // already checked
        if (neighbor.meta.currentPoint === p) continue;

        if (signedDistToPlane(p, neighbor.plane) > 0.0)
            getVisibleSet(p, neighbor, visible, horizon);
        else
            horizon.push(r);
    }
}

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

    generateFacetPlane(newFacet, points, centroid);

    return newFacet;
}

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

function createCentroid(points: Vector[], d: number): Vector
{
    // a point that will be internal from the very first simplex. Used to correctly orient new planes
    const centroid = points[0].slice();

    for (let j = 0; j < d; ++j) {
        for (let i = 1; i <= d; ++i) {
            centroid[j] += points[i][j];
        }
        centroid[j] /= d + 1;
    }

    return centroid;
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
    if (points.length === 0) return;
    const d = dim(points[0]);

    if (points.length <= d) {
        console.log(`A convex hull in ${d} dimensions requires at least ${d + 1} points.`);
    }

    const facets = createSimplex(points, d);
    for (let f of facets)
        f.meta = new FacetInfo();

    const centroid = createCentroid(points, d);

    // initial unprocessed point indices:
    const indices = [];
    for (let i = d + 1; i < points.length; ++i)
        indices.push(i);

    generateOutsideSets(indices, points, facets);

    // do not cache facets.length, as it will keep changing
    let done = false;

    // this extra loop should not be required
    while (!done) {
        done = true;
        for (let i = 0; i < facets.length; ++i) {
            const facet = facets[i];

            const p = getFurthestPoint(facet);
            if (p !== -1) {
                removeElementOutOfOrder(facet.meta.outsideSet, p);

                const V = [];
                const H = [];

                getVisibleSet(points[p], facet, V, H);

                const newFacets = connectHorizonRidges(points, p, H, centroid, d);

                for (let v of V) {
                    if (removeElementOutOfOrder(facets, v) <= i) --i;
                    generateOutsideSets(v.meta.outsideSet, points, newFacets);
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