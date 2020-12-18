import { Vector } from "../types";

export class Ridge
{
    verts: number[] = [];
    neighbor: Ridge;
    facet: Facet;

    constructor(facet: Facet)
    {
        this.facet = facet;
    }
}

export class Facet
{
    ridges: Ridge[] = [];
    plane: Vector;

    meta: any;
}