import { Facet, Ridge } from "./Geometry";
import { Vector } from "../types";
/**
 * Some shape construction code used internally.
 *
 * @author derschmale <http://www.derschmale.com>
 */
/**
 * Assign the neighbor from a set. Generally only used while constructing a new facet.
 *
 * @param facet The facet owning the ridge.
 * @param ridge The ridge for which to find the neighbor.
 * @param facets The set of facets to search.
 *
 * @ignore
 */
export declare function findNeighbor(facet: Facet, ridge: Ridge, facets: Facet[]): void;
/**
 * Generates a plane for a facet. Used when constructing new facets.
 * @param facet The facet to construct a new plane for.
 * @param points The general set of points the indices refer to.
 * @param centroid An optional point that's guaranteed to be behind the plane. This can be used to orient the face
 * if the vertex order is inconsistent.
 *
 * @ignore
 */
export declare function generateFacetPlane(facet: Facet, points: Vector[], dim: number, centroid?: Vector): void;
/**
 * Generates the ridges for a facet based on the vertices contained in it. If ridges are already present, they're
 * considered already built and only the missing ones will be added (used only when connecting ridges). In this
 * case, the order of vertices MUST match that of the present ridges.
 *
 * @ignore
 */
export declare function buildRidges(facet: Facet, facets: Facet[], dim: number): void;
