import { Vector } from "../types";
import { dim } from "../math/VecMath";
import { quickHull } from "../hull/QuickHull";
import { Facet } from "../geom/Geometry";

/**
 * 🎵🎶🎵  I lift you up to a higher dimension  🎶🎵🎶
 * 🎶🎵🎶     I'm gonna make you paraboloid     🎵🎶🎵
 *
 * @ignore
 */
function lift(points: Vector[], d: number): Vector[]
{
    const liftedDim = d + 1;
    const numPoints = points.length;
    // cache coherency
    const byteSize = liftedDim << 2; // every lifted point takes up 4 bytes per float
    const arr = new ArrayBuffer((numPoints + 1) * byteSize); // add room for one upper bound
    let i = 0;
    let bound = 0;
    let min = points[0].slice();
	let max = points[0].slice();

    const lifted = points.map(p =>
    {
        const f = new Float32Array(arr, i, liftedDim);
        let s = 0;

        for (let j = 0; j < d; ++j) {
            const e = p[j];
            // the random factor is cheating, but it seems to solve some precision errors if everything is on a grid
            f[j] = e;

            if (e < min[j])
				min[j] = e;

			if (e > max[j])
				max[j] = e;

            s += e * e;
        }

        f[d] = s;

        if (s > bound)
            bound = s;

        i += byteSize;

        return f;
    });

    // add a bounding point to increase robustness, will be filtered out on plane side test
    const boundPt = new Float32Array(arr, numPoints * byteSize, liftedDim);
    for (let i = 0; i < d; ++i) {
        boundPt[i] = (min[i] + max[i]) * .5;
    }
    boundPt[d] = bound;
    lifted.push(boundPt);

    return lifted;
}

/**
 * Calculates the Delaunay triangulation (or tetrahedralization) of a set of points.
 * @param points
 *
 * @author derschmale <http://www.derschmale.com>
 */
export function delaunay(points: Vector[]): Facet[]
{
    const d = dim(points[0]);
    const numPoints = points.length;

    if (numPoints === d + 1) {
        return quickHull(points);
    }

    const lifted = lift(points, d);
    const hull = quickHull(lifted);

    return hull.filter(f => {
        // remove all upwards facing faces and ridges
        if (f.plane[d] >= 0.0) {
            for (let r of f.ridges)
                if (r.neighbor) r.neighbor.neighbor = null;

            return false;
        }
        return true;
    });
}
