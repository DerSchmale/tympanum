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
export declare function dim(v: Vector): number;
/**
 * The dot (inner) product of two vectors.
 *
 * @ignore
 */
export declare function dot(v1: Vector, v2: Vector, dim?: number): number;
/**
 * Normalizes a vector.
 *
 * @ignore
 */
export declare function normalize(v: Vector): Vector;
/**
 * Normalizes a plane encoded as a vector as (normal, offset).
 *
 * @ignore
 */
export declare function normalizePlane(v: Vector): Vector;
/**
 * Creates a new N-dimensional matrix.
 * Indexing is [row][col].
 *
 * @ignore
 */
export declare function getSquareMatrix(dim: number): Vector[];
/**
 * Calculates the hyperplane that contains the given points. The amount of points defines the dimension of the
 * hyperplane.
 *
 * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
 * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
 *
 * @ignore
 */
export declare function hyperplaneFromPoints(p: Vector[], tgt?: Vector): Vector;
/**
 * Flips a vector.
 *
 * @ignore
 */
export declare function negate(v: Vector): Vector;
/**
 * Calculates the signed distance of a point to a plane.
 *
 * @ignore
 */
export declare function signedDistToPlane(point: Vector, plane: Vector, dim?: number): number;
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
export declare function intersectRayPlane(origin: Vector, dir: Vector, plane: Vector, dim: number, startsInside?: boolean): number;
/**
 * Calculates the inverse of a matrix
 * @ignore
 */
export declare function invertMatrix(mtx: Vector[], dim: number): Vector[];
/**
 * @ignore
 */
export declare function transformVector(mtx: Vector[], p: Vector, tgt: Vector, dim: number): Vector;
