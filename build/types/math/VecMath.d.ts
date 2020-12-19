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
export declare function signedDistToPlane(point: Vector, plane: Vector, dim: number): number;
