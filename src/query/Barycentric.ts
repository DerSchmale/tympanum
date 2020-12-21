import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
import { dim, getSquareMatrix, invertMatrix, transformVector } from "../math/VecMath";

/**
 * @ignore
 */
const mtxCache: Vector[][] = [];

/**
 * @ignore
 */
const tmpCache: Vector[] = [];

/**
 * Calculates the barycentric coordinates for a given point and a Facet. Every element of the coordinate is the
 * weight for the facet's vertex at the corresponding index.
 * @param position The position to calculate the barycentric coords for.
 * @param facet The Facet relative to which the coords are calculated.
 * @param points The points array indexed by Facet.
 * @param tgt An optional target to store the results. For dimension N, must be of length N+1.
 *
 * @author derschmale <http://www.derschmale.com>
 */
export function barycentricCoords(position: Vector, facet: Facet, points: Vector[], tgt?: Vector): Vector
{
    const d = dim(position);

    tgt = tgt || new Float32Array(d + 1);

    // we'll probably want to execute this function a lot of times, so let's make it efficient by not recreating these
    let mtx = mtxCache[d];
    if (!mtx) {
        mtxCache[d] = mtx = getSquareMatrix(d);
        tmpCache[d] = new Float32Array(d);
    }

    const tmp = tmpCache[d];
    const verts = facet.verts;
    const pN = points[verts[d]];

    // express all relative to pN
    for (let i = 0; i < d; ++i) {
        tmp[i] = position[i] - pN[i];

        const p = points[verts[i]];
        for (let j = 0; j < d; ++j) {
            // vectors are transposed!
            mtx[j][i] = p[j] - pN[j];
        }
    }

    invertMatrix(mtx, d);
    transformVector(mtx, tmp, tgt, d);
    tgt[d] = 1.0;

    for (let i = 0; i < d; ++i)
        tgt[d] -= tgt[i];

    return tgt;
}