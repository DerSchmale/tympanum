import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
/**
 * Calculates the Delaunay triangulation (or tetrahedralization) of a set of points.
 * @param points
 *
 * @author derschmale <http://www.derschmale.com>
 */
export declare function delaunay(points: Vector[]): Facet[];
