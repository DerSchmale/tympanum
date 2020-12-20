import { Vector } from "../types";
import { Facet } from "./Geometry";
/**
 * Creates an N-simplex from N+1 points. The dimension of the points are used to
 * define the dimension of the simplex.
 *
 * @param points An array of points. Only the first N+1 points are used.
 * @param indices An optional array of indices into points to define which points in the set are used.
 *
 * @author derschmale <http://www.derschmale.com>
 */
export declare function createSimplex(points: Vector[], indices?: number[]): Facet[];
