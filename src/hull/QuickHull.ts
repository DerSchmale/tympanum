import { Vector } from "../types";
import { dim, hyperplaneFromPoints, negate, signedDistToPlane } from "../math/VecMath";
import { removeElementOutOfOrder } from "../array/utils";
import { Facet, Ridge } from "../geom/Geometry";

class FacetInfo
{
    outsideSet: number[] = [];
    outsideDist: number[] = [];
    currentPoint: Vector;   // used to keep track of visibility tests
}

function createFacet(): Facet
{
    const f = new Facet();
    f.meta = new FacetInfo();
    return f;
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

function generatePlane(f: Facet, points: Vector[], centroid: Vector): void
{
    const verts = f.ridges.map(r => points[r.verts[0]]);

    let plane = f.plane = hyperplaneFromPoints(verts);
    // use an opposing point, which MUST be on the negative side of the plane

    if (signedDistToPlane(centroid, plane) > 0.0) {
        negate(plane);
        // flip ridges for consistency
        f.ridges.reverse();
        f.ridges.forEach(r => r.verts.reverse());
    }
}

function findNeighbor(facet: Facet, ridge: Ridge, facets: Facet[]): void
{
    const src = ridge.verts;
    const len = src.length;

    for (let f of facets) {
        for (let r of f.ridges) {
            // do not bother if we already found them
            if (r.neighbor) continue;

            let index = r.verts.indexOf(src[0]);
            let found = index >= 0;
            let i = 1;

            while (found && i < len) {
                index = (index + 1) % len;
                found = (r.verts[index]) == src[i];
                ++i;
            }

            if (found) {
                ridge.neighbor = r;
                r.neighbor = ridge;
                return;
            }
        }
    }
}

function createSimplex(points: Vector[], dim: number, centroid: Vector): Facet[]
{
    // construct facets from dim+1 points
    // 1D: [ [0], [1] ]
    // 2D: the set of 1D simplex facets for points 0, 1, 2
    //        [ [0, 1], [1, 2], [2, 0] ]
    // 3D: the set of 2D simplex facets for points 0, 1, 2, 3
    //    [
    //        [0, 1, 2]
    //        [1, 2, 3]
    //        [2, 3, 0]
    //        [3, 0, 2]
    //    ]
    // 4D: the set of 3D simplices facets for points 0, 1, 2, 3, 4
    //    [
    //        [0, 1, 2, 3]
    //        [1, 2, 3, 4]
    //        [2, 3, 4, 0]
    //        [3, 4, 0, 1]
    //        [4, 0, 1, 2]
    //    ]

    // TODO: We can find a better initial data set, similar to QHull:
    //  find the minX, maxX points, and iteratively extend with furthest point
    //  (starts with signed dist to line, then to plane in 3D, then to hyperplane in 4D)
    //  This is the same sort of logic of the base Quickhull algorithm, so maybe it's not that much of an improvement

    const facets: Facet[] = [];
    // the facet is made up of dim + 1 points, so cycle through these
    const numVerts = dim + 1;

    // Dimension of simplex = dim
    // Number of facets = dim + 1
    //  2D: 3 lines to a triangle
    //  3D: 4 triangles to a tetrahedron
    // Number of ridges per facet = dim
    //  2D: 2 vertices to a segment
    //  3D: 3 edges to a triangle
    // Number of vertices per ridge = dim - 1
    //  2D: 1 point per vertex
    //  3D: 2 vertices to an edge

    let verts = [];

    for (let i = 0; i <= dim; ++i) {
        let f = createFacet();

        // collect all verts for this facet
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

        // we need dim + 1 points
        generatePlane(f, points, centroid);

        facets.push(f);
    }

    return facets;
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
    const newFacet = createFacet();

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

    generatePlane(newFacet, points, centroid);

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

    // a point that will be internal from the very first simplex. Used to correctly orient new planes
    const centroid = points[0].slice();

    for (let j = 0; j < d; ++j) {
        for (let i = 1; i <= d; ++i) {
            centroid[j] += points[i][j];
        }
        centroid[j] /= d + 1;
    }

    const facets = createSimplex(points, d, centroid);

    // initial unprocessed point indices:
    const indices = [];

    for (let i = d + 1; i < points.length; ++i) {
        indices.push(i);
    }

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

    facets.forEach(f => f.meta = null);

    return facets;
}