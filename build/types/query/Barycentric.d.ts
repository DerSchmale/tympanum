import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
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
export declare function barycentricCoords(position: Vector, facet: Facet, points: Vector[], tgt?: Vector): Vector;
