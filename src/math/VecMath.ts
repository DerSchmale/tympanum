import { Vector } from "../types";

/**
 * Some basic N-dimensional vector math.
 * @author derschmale <http://www.derschmale.com>
 */

/**
 * Returns the dimension of a Vector.
 *
 * @ignore
 */
export function dim(v: Vector): number
{
    return v.length;
}

/**
 * The dot (inner) product of two vectors.
 *
 * @ignore
 */
export function dot(v1: Vector, v2: Vector, dim: number = -1): number
{
    if (dim === -1)
        dim = v1.length;

    let d = v1[0] * v2[0];

    for (let i = 1; i < dim; ++i)
        d += v1[i] * v2[i];

    return d;
}

/**
 * Normalizes a vector.
 *
 * @ignore
 */
export function normalize(v: Vector, dim?: number): Vector
{
    dim = dim ?? v.length;
    let d = v[0] * v[0];

    for (let i = 1; i < dim; ++i)
        d += v[i] * v[i];

    if (d === 0.0) return v;
    d = 1.0 / Math.sqrt(d);

    for (let i = 0; i < dim; ++i)
        v[i] *= d;

    return v;
}

/**
 * Normalizes a plane encoded as a vector as (normal, offset).
 *
 * @ignore
 */
export function normalizePlane(v: Vector): Vector
{
    const len = v.length;
    let d = v[0] * v[0];

    for (let i = 1; i < len - 1; ++i)
        d += v[i] * v[i];

    if (d === 0.0) return v;
    d = 1.0 / Math.sqrt(d);

    for (let i = 0; i < len; ++i)
        v[i] *= d;

    return v;
}

/**
 * Generates the cofactor matrix.
 *
 * @ignore
 */
function cofactor(mat: Vector[], tgt: Vector[], row: number, col: number, dim: number)
{
    let i = 0, j = 0;

    for (let r = 0; r < dim; ++r) {
        if (r === row) continue;

        j = 0;

        for (let c = 0; c < dim; ++c) {
            if (c === col) continue;

            tgt[i][j] = mat[r][c];
            ++j;
        }

        ++i;
    }
}

/**
 * Creates a new N-dimensional matrix.
 * Indexing is [row][col].
 *
 * @ignore
 */
export function getSquareMatrix(dim: number): Vector[]
{
    const data = new ArrayBuffer(dim * dim * 4);
    const mtx = [];

    for (let i = 0; i < dim; ++i)
        mtx[i] = new Float32Array(data, i * dim << 2, dim);

    return mtx;
}

/**
 * Calculates the determinant for a matrix.
 *
 * @ignore
 */
function det(v: Vector[], dim: number): number
{
    if (dim === 1) {
        return v[0][0];
    }
    // just early base cases for efficiency, these wouldn't really be necessary
    else if (dim === 2) {
        return v[0][0] * v[1][1] - v[0][1] * v[1][0];
    }
    else if (dim === 3) {
        return v[0][0] * v[1][1] * v[2][2] + v[0][1] * v[1][2] * v[2][0] + v[0][2] * v[1][0] * v[2][1]
            - v[0][2] * v[1][1] * v[2][0] - v[0][1] * v[1][0] * v[2][2] - v[0][0] * v[1][2] * v[2][1];
    }
    else {
        let s = 1;
        let d = 0;
        let sub = getSquareMatrix(dim - 1);

        for (let i = 0; i < dim; ++i) {
            cofactor(v, sub, 0, i, dim);
            d += s * v[0][i] * det(sub, dim - 1);
            s = -s;
        }
        return d;
    }
}

/**
 * Calculates the generalized cross product.
 *
 * @ignore
 */
function generalizedCross(v: Vector[], tgt?: Vector): Vector
{
    const dim = v[0].length;

    tgt = tgt || new Float32Array(dim);

    if (dim === 3) {
        const ar = v[0], br = v[1];
        const ax = ar[0], ay = ar[1], az = ar[2];
        const bx = br[0], by = br[1], bz = br[2];

        tgt[0] = ay * bz - az * by;
        tgt[1] = az * bx - ax * bz;
        tgt[2] = ax * by - ay * bx;
    }
    else {
        // so the generalized cross product is the determinant of:
        // | v00 v01 v02 ... |
        // | v10 v11 v12 ... |
        // | v20 v21 v22 ... |
        // | ... ... ... ... |
        // | e0  e1  e2  ... |
        // with eN being the n-th standard basis vector (1, 0, 0), (0, 1, 0)
        // IE: the eN elements are basically "selectors" for each target vector's element

        let sign = dim % 2 ? -1 : 1;
        let sub = getSquareMatrix(dim - 1);

        for (let i = 0; i < dim; ++i) {
            // use the last row because those would contain the basis vectors
            // pass in v as a matrix, which will look like the actual matrix
            cofactor(v, sub, dim - 1, i, dim);
            tgt[i] = sign * det(sub, dim - 1);
            sign = -sign;
        }
    }

    return tgt;
}

/**
 * Calculates the hyperplane that contains the given points. The amount of points defines the dimension of the
 * hyperplane.
 *
 * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
 * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
 *
 * @ignore
 */
export function hyperplaneFromPoints(p: Vector[], tgt?: Vector)
{
    const dim = p.length;
    const v0 = p[0];
    const vecs = [];

    tgt = tgt || new Float32Array(dim + 1);

    for (let i = 1; i < dim; ++i) {
        const pt = p[i];
        const v = [];

        for (let j = 0; j < dim; ++j)
            v[j] = v0[j] - pt[j];

        vecs.push(v);
    }

    // calculate normal for hyperplane
    generalizedCross(vecs, tgt);

    tgt[dim] = -dot(v0, tgt, dim);
    normalizePlane(tgt);

    return tgt;
}

/**
 * Flips a vector.
 *
 * @ignore
 */
export function negate(v: Vector): Vector
{
    const dim = v.length;
    for (let i = 0; i < dim; ++i)
        v[i] = -v[i];
    return v;
}

/**
 * Calculates the signed distance of a point to a plane.
 *
 * @ignore
 */
export function signedDistToPlane(point: Vector, plane: Vector, dim: number = -1)
{
    if (dim === -1)
        dim = point.length;

    let d = plane[dim];

    for (let i = 0; i < dim; ++i)
        d += point[i] * plane[i];

    return d;
}

/**
 * Calculates the intersection with a ridge's hyperplane and a ray
 *
 * @param origin The origin of the ray.
 * @param dir The direction of the ray. Doesn't need to be normalized. When testing segments, this is (end - start).
 * @param plane The plane to test against.
 * @param dim The dimension.
 * @param startsInside Set to true if we're testing for intersections of planes of a convex solid and the start
 * point is inside. Used for early rejection tests if the planes are facing away from the ray.
 *
 * @ignore
 */
export function intersectRayPlane(origin: Vector, dir: Vector, plane: Vector, dim: number, startsInside: boolean = false): number
{
    // assuming vectors are all normalized
    const denom = dot(plane, dir, dim);
    const abs = Math.abs(denom);

    // not parallel + must be traveling in the same direction if it starts inside
    if (abs > 0.0 && (!startsInside || denom > 0.0)) {
        return -signedDistToPlane(origin, plane, dim) / denom;
    }

    return -1;
}

/**
 * Calculates the adjoint of a matrix
 *
 * @ignore
 */
function adjointMatrix(mtx: Vector[], tgt: Vector[], dim: number): Vector[]
{
    if (dim === 1)
        return [ [ 1 ] ];

    // temp is used to store cofactors of A[][]
    let sign = 1;
    let cof = getSquareMatrix(dim);

    for (let i = 0; i < dim; ++i) {
        for (let j = 0; j < dim; ++j) {
            cofactor(mtx, cof, i, j, dim);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            tgt[j][i] = (sign) * (det(cof, dim - 1));
        }
    }

    return tgt;
}

/**
 * Calculates the inverse of a matrix
 * @ignore
 */
export function invertMatrix(mtx: Vector[], dim: number): Vector[]
{
    // custom cases for efficiency
    if (dim === 2) {
        const m00 = mtx[0][0], m01 = mtx[0][1];
        const m10 = mtx[1][0], m11 = mtx[1][1];
        const determinant = m00 * m11 - m01 * m10;
        if (determinant === 0.0) return null;
        const rcpDet = 1.0 / determinant;
        mtx[0][0] = m11 * rcpDet;
        mtx[0][1] = -m01 * rcpDet;
        mtx[1][0] = -m10 * rcpDet;
        mtx[1][1] = m00 * rcpDet;
    }
    else if (dim === 3) {
        const m00 = mtx[0][0], m01 = mtx[0][1], m02 = mtx[0][2];
        const m10 = mtx[1][0], m11 = mtx[1][1], m12 = mtx[1][2];
        const m20 = mtx[2][0], m21 = mtx[2][1], m22 = mtx[2][2];

        const determinant = m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
        if (determinant === 0.0) return null;
        const rcpDet = 1.0 / determinant;

        mtx[0][0] = (m11 * m22 - m12 * m21) * rcpDet;
        mtx[0][1] = (m02 * m21 - m01 * m22) * rcpDet;
        mtx[0][2] = (m01 * m12 - m02 * m11) * rcpDet;
        mtx[1][0] = (m12 * m20 - m10 * m22) * rcpDet;
        mtx[1][1] = (m00 * m22 - m02 * m20) * rcpDet;
        mtx[1][2] = (m02 * m10 - m00 * m12) * rcpDet;
        mtx[2][0] = (m10 * m21 - m11 * m20) * rcpDet;
        mtx[2][1] = (m01 * m20 - m00 * m21) * rcpDet;
        mtx[2][2] = (m00 * m11 - m01 * m10) * rcpDet;
        return this;
    }
    else {
        // There are faster ways of doing this, but for now, it'll do
        let determinant = det(mtx, dim);
        if (determinant === 0) {
            return null;
        }

        const rcpDet = 1.0 / determinant;

        // Find adjoint
        let adj = adjointMatrix(mtx, getSquareMatrix(dim), dim);

        // inverse(A) = adjoint/determinant
        for (let i = 0; i < dim; ++i)
            for (let j = 0; j < dim; j++)
                mtx[i][j] = adj[i][j] * rcpDet;

    }
    return mtx;
}

/**
 * @ignore
 */
export function transformVector(mtx: Vector[], p: Vector, tgt: Vector, dim: number)
{
    for (let i = 0; i < dim; ++i) {
        tgt[i] = 0;
        for (let j = 0; j < dim; ++j) {
            tgt[i] += p[j] * mtx[i][j];
        }
    }

    return tgt;
}