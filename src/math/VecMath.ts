import { Vector } from "../types";

export function dim(v: Vector): number
{
    return v.length;
}

export function dot(v1: Vector, v2: Vector): number
{
    const dim = v1.length;
    let d = v1[0] * v2[0];

    for (let i = 1; i < dim; ++i)
        d += v1[i] * v2[i];

    return d;
}

export function normalize(v: Vector): Vector
{
    const dim = v.length;
    let d = v[0] * v[0];

    for (let i = 1; i < dim; ++i)
        d += v[i] * v[i];

    if (d === 0.0) return v;
    d = 1.0 / Math.sqrt(d);

    for (let i = 0; i < dim; ++i)
        v[i] *= d;

    return v;
}

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

function cofactor(mat: Vector[], tgt: Vector[], p: number, q: number, dim: number)
{
    let i = 0, j = 0;

    for (let row = 0; row < dim; row++) {
        if (row === p) continue;

        for (let col = 0; col < dim; col++) {
            if (col === q) continue;

            tgt[i][j] = mat[row][col];
            ++j;
        }

        ++i;
    }
}

function getSquareMatrix(dim: number): Vector[]
{
    let sub = [];

    for (let i = 0; i < dim; ++i)
        sub[i] = new Float32Array(dim);

    return sub;
}

// if this is a determinant for a submatrix, ij are the indices of the parent
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
        let sub = getSquareMatrix(dim);

        for (let i = 0; i < dim; ++i) {
            cofactor(v, sub, 0, i, dim);
            d += s * v[0][i] * det(sub, dim - 1);
            s = -s;
        }
        return d;
    }
}

// there's probably waaaay faster algos to do this
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

        let sign = dim % 2? -1 : 1;
        let sub = getSquareMatrix(dim);

        for (let i = 0; i < dim; ++i) {
            cofactor(v, sub, dim - 1, i, dim);
            tgt[i] = sign * det(sub, dim - 1);
            sign = -sign;
        }
    }

    return tgt;
}

export function hyperplaneFromPoints(p: Vector[], tgt?: Vector) {
    const dim = p.length;
    const v0 = p[0];
    const vecs = [];

    tgt = tgt || new Float32Array(dim + 1);

    for (let i = 1; i < dim; ++i) {
        const pt = p[i];
        const v = [];

        for (let j = 0; j < dim; ++j) {
            v[j] = pt[j] - v0[j];
        }

        vecs.push(v);
    }

    // calculate normal for hyperplane
    generalizedCross(vecs, tgt);
    // calculate offset
    tgt[dim] = -dot(v0, tgt);
    // not sure if this is necessary
    normalizePlane(tgt);
    return tgt;
}

export function negate(v: Vector): Vector
{
    const dim = v.length;
    for (let i = 0; i < dim; ++i)
        v[i] = -v[i];
    return v;
}
export function signedDistToPlane(v: Vector, p: Vector)
{
    const dim = v.length;
    let d = p[dim];

    for (let i = 0; i < dim; ++i)
        d += v[i] * p[i];

    return d;
}