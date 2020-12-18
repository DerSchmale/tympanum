var TYMP = (function (exports) {
    'use strict';

    function dim(v) {
        return v.length;
    }
    function dot(v1, v2) {
        var dim = v1.length;
        var d = v1[0] * v2[0];
        for (var i = 1; i < dim; ++i)
            d += v1[i] * v2[i];
        return d;
    }
    function normalizePlane(v) {
        var len = v.length;
        var d = v[0] * v[0];
        for (var i = 1; i < len - 1; ++i)
            d += v[i] * v[i];
        if (d === 0.0)
            return v;
        d = 1.0 / Math.sqrt(d);
        for (var i = 0; i < len; ++i)
            v[i] *= d;
        return v;
    }
    function cofactor(mat, tgt, p, q, dim) {
        var i = 0, j = 0;
        for (var row = 0; row < dim; row++) {
            if (row === p)
                continue;
            for (var col = 0; col < dim; col++) {
                if (col === q)
                    continue;
                tgt[i][j] = mat[row][col];
                ++j;
            }
            ++i;
        }
    }
    function getSquareMatrix(dim) {
        var sub = [];
        for (var i = 0; i < dim; ++i)
            sub[i] = new Float32Array(dim);
        return sub;
    }
    // if this is a determinant for a submatrix, ij are the indices of the parent
    function det(v, dim) {
        if (dim === 1) {
            return v[0][0];
        }
        // just early base cases for efficiency, these wouldn't really be necessary
        else if (dim === 2) {
            return v[0][0] * v[1][1] - v[0][1] * v[1][0];
        }
        else if (dim === 3) {
            return v[0][0] * v[1][1] * v[2][2] + v[0][1] * v[1][2] * v[2][0] + v[0][2] * v[1][0] * v[2][1]
                - v[0][2] * v[1][1] * v[2][0] - v[0][1] * v[1][0] * v[2][2] - v[0][0] * v[1][2] * v[2][1];
        }
        else {
            var s = 1;
            var d = 0;
            var sub = getSquareMatrix(dim);
            for (var i = 0; i < dim; ++i) {
                cofactor(v, sub, 0, i, dim);
                d += s * v[0][i] * det(sub, dim - 1);
                s = -s;
            }
            return d;
        }
    }
    // there's probably waaaay faster algos to do this
    function generalizedCross(v, tgt) {
        var dim = v[0].length;
        tgt = tgt || new Float32Array(dim);
        if (dim === 3) {
            var ar = v[0], br = v[1];
            var ax = ar[0], ay = ar[1], az = ar[2];
            var bx = br[0], by = br[1], bz = br[2];
            tgt[0] = ay * bz - az * by;
            tgt[1] = az * bx - ax * bz;
            tgt[2] = ax * by - ay * bx;
        }
        else {
            // so the generalized cross product is the determinant of:
            // | v00 v01 v02 ... |
            // | v10 v11 v12 ... |
            // | v20 v21 v22 ... |
            // | ... ... ... ... |
            // | e0  e1  e2  ... |
            // with eN being the n-th standard basis vector (1, 0, 0), (0, 1, 0)
            // IE: the eN elements are basically "selectors" for each target vector's element
            var sign = dim % 2 ? -1 : 1;
            var sub = getSquareMatrix(dim);
            for (var i = 0; i < dim; ++i) {
                cofactor(v, sub, dim - 1, i, dim);
                tgt[i] = sign * det(sub, dim - 1);
                sign = -sign;
            }
        }
        return tgt;
    }
    function hyperplaneFromPoints(p, tgt) {
        var dim = p.length;
        var v0 = p[0];
        var vecs = [];
        tgt = tgt || new Float32Array(dim + 1);
        for (var i = 1; i < dim; ++i) {
            var pt = p[i];
            var v = [];
            for (var j = 0; j < dim; ++j) {
                v[j] = pt[j] - v0[j];
            }
            vecs.push(v);
        }
        // calculate normal for hyperplane
        generalizedCross(vecs, tgt);
        // calculate offset
        tgt[dim] = -dot(v0, tgt);
        // not sure if this is necessary
        normalizePlane(tgt);
        return tgt;
    }
    function negate(v) {
        var dim = v.length;
        for (var i = 0; i < dim; ++i)
            v[i] = -v[i];
        return v;
    }
    function signedDistToPlane(v, p) {
        var dim = v.length;
        var d = p[dim];
        for (var i = 0; i < dim; ++i)
            d += v[i] * p[i];
        return d;
    }

    function removeElementOutOfOrder(target, elm) {
        var last = target.pop();
        if (last === elm) {
            return target.length;
        }
        else {
            var index = target.indexOf(elm);
            if (index === -1) {
                target.push(last);
                throw new Error("Removing component that's not present");
            }
            target[index] = last;
            return index;
        }
    }

    var Ridge = /** @class */ (function () {
        function Ridge(facet) {
            this.verts = [];
            this.facet = facet;
        }
        return Ridge;
    }());
    var Facet = /** @class */ (function () {
        function Facet() {
            this.ridges = [];
        }
        return Facet;
    }());
    var FacetInfo = /** @class */ (function () {
        function FacetInfo() {
            this.outsideSet = [];
            this.outsideDist = [];
        }
        return FacetInfo;
    }());
    function createFacet() {
        var f = new Facet();
        f.meta = new FacetInfo();
        return f;
    }
    function getFurthestPoint(facet) {
        var _a = facet.meta, outsideSet = _a.outsideSet, outsideDist = _a.outsideDist;
        var len = outsideSet.length;
        if (len === 0)
            return -1;
        var p = outsideSet[0];
        var dist = outsideDist[0];
        for (var i = 1; i < len; ++i) {
            if (outsideDist[i] > dist) {
                dist = outsideDist[i];
                p = outsideSet[i];
            }
        }
        return p;
    }
    function generatePlane(f, points, centroid) {
        var verts = f.ridges.map(function (r) { return points[r.verts[0]]; });
        var plane = f.plane = hyperplaneFromPoints(verts);
        // use an opposing point, which MUST be on the negative side of the plane
        if (signedDistToPlane(centroid, plane) > 0.0) {
            negate(plane);
            // flip ridges for consistency
            f.ridges.reverse();
            f.ridges.forEach(function (r) { return r.verts.reverse(); });
        }
    }
    function findNeighbor(facet, ridge, facets) {
        var src = ridge.verts;
        var len = src.length;
        for (var _i = 0, facets_1 = facets; _i < facets_1.length; _i++) {
            var f = facets_1[_i];
            for (var _a = 0, _b = f.ridges; _a < _b.length; _a++) {
                var r = _b[_a];
                // do not bother if we already found them
                if (r.neighbor)
                    continue;
                var index = r.verts.indexOf(src[0]);
                var found = index >= 0;
                var i = 1;
                while (found && i < len) {
                    index = (index + 1) % len;
                    found = (r.verts[index]) == src[i];
                    ++i;
                }
                if (found) {
                    ridge.neighbor = r;
                    r.neighbor = ridge;
                    return;
                }
            }
        }
    }
    function createSimplex(points, dim, centroid) {
        // construct facets from dim+1 points
        // 1D: [ [0], [1] ]
        // 2D: the set of 1D simplex facets for points 0, 1, 2
        //        [ [0, 1], [1, 2], [2, 0] ]
        // 3D: the set of 2D simplex facets for points 0, 1, 2, 3
        //    [
        //        [0, 1, 2]
        //        [1, 2, 3]
        //        [2, 3, 0]
        //        [3, 0, 2]
        //    ]
        // 4D: the set of 3D simplices facets for points 0, 1, 2, 3, 4
        //    [
        //        [0, 1, 2, 3]
        //        [1, 2, 3, 4]
        //        [2, 3, 4, 0]
        //        [3, 4, 0, 1]
        //        [4, 0, 1, 2]
        //    ]
        // TODO: We can find a better initial data set, similar to QHull:
        //  find the minX, maxX points, and iteratively extend with furthest point
        //  (starts with signed dist to line, then to plane in 3D, then to hyperplane in 4D)
        //  This is the same sort of logic of the base Quickhull algorithm, so maybe it's not that much of an improvement
        var facets = [];
        // the facet is made up of dim + 1 points, so cycle through these
        var numVerts = dim + 1;
        // Dimension of simplex = dim
        // Number of facets = dim + 1
        //  2D: 3 lines to a triangle
        //  3D: 4 triangles to a tetrahedron
        // Number of ridges per facet = dim
        //  2D: 2 vertices to a segment
        //  3D: 3 edges to a triangle
        // Number of vertices per ridge = dim - 1
        //  2D: 1 point per vertex
        //  3D: 2 vertices to an edge
        var verts = [];
        for (var i = 0; i <= dim; ++i) {
            var f = createFacet();
            // collect all verts for this facet
            for (var v = 0; v < dim; ++v) {
                verts[v] = (i + v) % numVerts;
            }
            for (var r = 0; r < dim; ++r) {
                var ridge = f.ridges[r] = new Ridge(f);
                for (var v = 0; v < dim - 1; ++v) {
                    ridge.verts[v] = verts[(r + v) % dim];
                }
                // find neighbours from already generated sets
                // there's probably an analytical way to do this
                findNeighbor(f, ridge, facets);
            }
            // we need dim + 1 points
            generatePlane(f, points, centroid);
            facets.push(f);
        }
        return facets;
    }
    // [face][setIndex][0/1] : 0 = index, 1 = signed distance to facet plane
    // assign points to the outside set of a collection of faces
    function generateOutsideSets(indices, points, facets) {
        var outsideSets = facets.map(function (_) { return []; });
        var len = indices.length;
        for (var i = 0; i < len; ++i) {
            var index = indices[i];
            var p = points[index];
            for (var _i = 0, facets_2 = facets; _i < facets_2.length; _i++) {
                var f = facets_2[_i];
                var dist = signedDistToPlane(p, f.plane);
                if (dist > 0) {
                    var meta = f.meta;
                    meta.outsideSet.push(index);
                    meta.outsideDist.push(dist);
                    break;
                }
            }
        }
        // any point now not in an outside set is inside the hull
        return outsideSets;
    }
    function getVisibleSet(p, facet, visible, horizon) {
        visible.push(facet);
        facet.meta.currentPoint = p;
        for (var _i = 0, _a = facet.ridges; _i < _a.length; _i++) {
            var r = _a[_i];
            var neighbor = r.neighbor.facet;
            // already checked
            if (neighbor.meta.currentPoint === p)
                continue;
            if (signedDistToPlane(p, neighbor.plane) > 0.0)
                getVisibleSet(p, neighbor, visible, horizon);
            else
                horizon.push(r);
        }
    }
    function attachNewFacet(ridge, p, points, facets, centroid, dim) {
        // in 2D, we simply need to create 1 new facet (line) from old ridge to p
        var newFacet = createFacet();
        // collect all verts for this facet, which is the horizon ridge + this point
        var verts = ridge.verts.slice();
        verts.push(p);
        // the horizon ridge is part of the new facet, and gets to keep its neighbor
        newFacet.ridges.push(ridge);
        ridge.facet = newFacet;
        // dim + 1 ridges (3 edges to a triangle in 2D, 4 faces to a tetrahedron in 3D)
        // start with 1, since 0 would be the same as the already existing ridge above
        for (var r = 1; r < dim; ++r) {
            var ridge_1 = new Ridge(newFacet);
            for (var v = 0; v < dim - 1; ++v)
                ridge_1.verts[v] = verts[(r + v) % dim];
            // so we only need to search for neighbours in the newly generated facets, only the horizons are attached to
            // the old ones
            findNeighbor(newFacet, ridge_1, facets);
            newFacet.ridges.push(ridge_1);
        }
        generatePlane(newFacet, points, centroid);
        return newFacet;
    }
    function connectHorizonRidges(points, index, H, centroid, dim) {
        var newFacets = [];
        // link horizon ridges with new point
        for (var _i = 0, H_1 = H; _i < H_1.length; _i++) {
            var ridge = H_1[_i];
            var newFacet = attachNewFacet(ridge, index, points, newFacets, centroid, dim);
            newFacets.push(newFacet);
        }
        return newFacets;
    }
    /**
     * QuickHull implements the algorithm of the same name, based on the original paper by Barber, Dobkin and Huhdanpaa.
     * We're not interested in 0- or 1-dimensional cases (the latter can simply be sorted points). QuickHull returns a
     * set of indices into the original point list so we can map it to a different original ata set (fe: points may be a
     * mapping for position vectors on some scene graph object).
     *
     * @author derschmale <http://www.derschmale.com>
     */
    function quickHull(points) {
        if (points.length === 0)
            return;
        var d = dim(points[0]);
        if (points.length <= d) {
            console.log("A convex hull in " + d + " dimensions requires at least " + (d + 1) + " points.");
        }
        // a point that will be internal from the very first simplex. Used to correctly orient new planes
        var centroid = points[0].slice();
        for (var j = 0; j < d; ++j) {
            for (var i = 1; i <= d; ++i) {
                centroid[j] += points[i][j];
            }
            centroid[j] /= d + 1;
        }
        var facets = createSimplex(points, d, centroid);
        // initial unprocessed point indices:
        var indices = [];
        for (var i = d + 1; i < points.length; ++i) {
            indices.push(i);
        }
        generateOutsideSets(indices, points, facets);
        // do not cache facets.length, as it will keep changing
        var done = false;
        // this extra loop should not be required
        while (!done) {
            done = true;
            for (var i = 0; i < facets.length; ++i) {
                var facet = facets[i];
                var p = getFurthestPoint(facet);
                if (p !== -1) {
                    removeElementOutOfOrder(facet.meta.outsideSet, p);
                    var V = [];
                    var H = [];
                    getVisibleSet(points[p], facet, V, H);
                    var newFacets = connectHorizonRidges(points, p, H, centroid, d);
                    for (var _i = 0, V_1 = V; _i < V_1.length; _i++) {
                        var v = V_1[_i];
                        if (removeElementOutOfOrder(facets, v) <= i)
                            --i;
                        generateOutsideSets(v.meta.outsideSet, points, newFacets);
                        if (v.meta.outsideSet.length > 0)
                            done = false;
                    }
                    facets.push.apply(facets, newFacets);
                }
            }
        }
        facets.forEach(function (f) { return f.meta = null; });
        return facets;
    }

    exports.quickHull = quickHull;

    Object.defineProperty(exports, '__esModule', { value: true });

    return exports;

}({}));
