import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
/**
 * QuickHull implements the algorithm of the same name, based on the original paper by Barber, Dobkin and Huhdanpaa.
 * We're not interested in 0- or 1-dimensional cases (the latter can simply be the extent of the point values).
 * QuickHull returns a set of indices into the original point list so we can map it to a different original ata set
 * (fe: points may be a mapping for position vectors on some scene graph object).
 *
 * @author derschmale <http://www.derschmale.com>
 */
export declare function quickHull(points: Vector[]): Facet[];
