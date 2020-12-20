import { Vector } from "../types";
import { Facet } from "../geom/Geometry";
/**
 * Performs the visibility walk algorithm to find the Facet containing the given position. This should only be used
 * on Delaunay triangulations, as other triangulations are not guaranteed to resolve to a solution.
 * @param position The position to search for.
 * @param facets The facets to search
 * @param points The points indexed by the facets.
 * @param startFacet An optional facet to start the search. If none is provided, an initial search
 * estimate is made, but this is not guaranteed to be a performance improvement.
 */
export declare function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], startFacet?: Facet): Facet;
/**
 * Performs the visibility walk algorithm to find the Facet containing the given position. This should only be used
 * on Delaunay triangulations, as other triangulations are not guaranteed to resolve to a solution.
 * @param position The position to search for.
 * @param facets The facets to search
 * @param points The points indexed by the facets.
 * @param estimate If true, searches the facets to find an initial match.
 */
export declare function visibilityWalk(position: Vector, facets: Facet[], points: Vector[], estimate?: boolean): Facet;
