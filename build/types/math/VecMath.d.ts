import { Vector } from "../types";
export declare function dim(v: Vector): number;
export declare function dot(v1: Vector, v2: Vector): number;
export declare function normalize(v: Vector): Vector;
export declare function normalizePlane(v: Vector): Vector;
export declare function hyperplaneFromPoints(p: Vector[], tgt?: Vector): Vector;
export declare function negate(v: Vector): Vector;
export declare function signedDistToPlane(v: Vector, p: Vector): number;
