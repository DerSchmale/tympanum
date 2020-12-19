import { Vector } from "../types";
import { Facet } from "./Geometry";
/**
 * Creates an N-simplex from N+1 points.
 *
 * @param points An array of points. Only the first N+1 points are used.
 * @param dim The dimension of the simplex.
 * @param indices An optional array of indices into points to remap which points are used
 *
 * @author derschmale <http://www.derschmale.com>
 */
export declare function createSimplex(points: Vector[], dim: number, indices?: number[]): Facet[];
