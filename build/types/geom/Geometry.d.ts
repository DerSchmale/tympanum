import { Vector } from "../types";
export declare class Ridge {
    verts: number[];
    neighbor: Ridge;
    facet: Facet;
    constructor(facet: Facet);
}
export declare class Facet {
    ridges: Ridge[];
    plane: Vector;
    meta: any;
}
