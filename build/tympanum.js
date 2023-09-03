var TYMP = (function (exports) {
    'use strict';

    /**
     * Some basic N-dimensional vector math.
     * @author derschmale <http://www.derschmale.com>
     */
    /**
     * Returns the dimension of a Vector.
     *
     * @ignore
     */
    function dim(v) {
        return v.length;
    }
    /**
     * The dot (inner) product of two vectors.
     *
     * @ignore
     */
    function dot(v1, v2, dim) {
        if (dim === void 0) { dim = -1; }
        if (dim === -1)
            dim = v1.length;
        var d = v1[0] * v2[0];
        for (var i = 1; i < dim; ++i)
            d += v1[i] * v2[i];
        return d;
    }
    /**
     * Normalizes a plane encoded as a vector as (normal, offset).
     *
     * @ignore
     */
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
    /**
     * Generates the cofactor matrix.
     *
     * @ignore
     */
    function cofactor(mat, tgt, row, col, dim) {
        var i = 0, j = 0;
        for (var r = 0; r < dim; ++r) {
            if (r === row)
                continue;
            j = 0;
            for (var c = 0; c < dim; ++c) {
                if (c === col)
                    continue;
                tgt[i][j] = mat[r][c];
                ++j;
            }
            ++i;
        }
    }
    /**
     * Creates a new N-dimensional matrix.
     * Indexing is [row][col].
     *
     * @ignore
     */
    function getSquareMatrix(dim) {
        var data = new ArrayBuffer(dim * dim * 4);
        var mtx = [];
        for (var i = 0; i < dim; ++i)
            mtx[i] = new Float32Array(data, i * dim << 2, dim);
        return mtx;
    }
    /**
     * Calculates the determinant for a matrix.
     *
     * @ignore
     */
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
            var sub = getSquareMatrix(dim - 1);
            for (var i = 0; i < dim; ++i) {
                cofactor(v, sub, 0, i, dim);
                d += s * v[0][i] * det(sub, dim - 1);
                s = -s;
            }
            return d;
        }
    }
    /**
     * Calculates the generalized cross product.
     *
     * @ignore
     */
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
            var sub = getSquareMatrix(dim - 1);
            for (var i = 0; i < dim; ++i) {
                // use the last row because those would contain the basis vectors
                // pass in v as a matrix, which will look like the actual matrix
                cofactor(v, sub, dim - 1, i, dim);
                tgt[i] = sign * det(sub, dim - 1);
                sign = -sign;
            }
        }
        return tgt;
    }
    /**
     * Calculates the hyperplane that contains the given points. The amount of points defines the dimension of the
     * hyperplane.
     *
     * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
     * 🎵🎶🎵  Hypah Hypah!  🎶🎵🎶
     *
     * @ignore
     */
    function hyperplaneFromPoints(p, tgt) {
        var dim = p.length;
        var v0 = p[0];
        var vecs = [];
        tgt = tgt || new Float32Array(dim + 1);
        for (var i = 1; i < dim; ++i) {
            var pt = p[i];
            var v = [];
            for (var j = 0; j < dim; ++j)
                v[j] = v0[j] - pt[j];
            vecs.push(v);
        }
        // calculate normal for hyperplane
        generalizedCross(vecs, tgt);
        tgt[dim] = -dot(v0, tgt, dim);
        normalizePlane(tgt);
        return tgt;
    }
    /**
     * Flips a vector.
     *
     * @ignore
     */
    function negate(v) {
        var dim = v.length;
        for (var i = 0; i < dim; ++i)
            v[i] = -v[i];
        return v;
    }
    /**
     * Calculates the signed distance of a point to a plane.
     *
     * @ignore
     */
    function signedDistToPlane(point, plane, dim) {
        if (dim === void 0) { dim = -1; }
        if (dim === -1)
            dim = point.length;
        var d = plane[dim];
        for (var i = 0; i < dim; ++i)
            d += point[i] * plane[i];
        return d;
    }
    /**
     * Calculates the intersection with a ridge's hyperplane and a ray
     *
     * @param origin The origin of the ray.
     * @param dir The direction of the ray. Doesn't need to be normalized. When testing segments, this is (end - start).
     * @param plane The plane to test against.
     * @param dim The dimension.
     * @param startsInside Set to true if we're testing for intersections of planes of a convex solid and the start
     * point is inside. Used for early rejection tests if the planes are facing away from the ray.
     *
     * @ignore
     */
    function intersectRayPlane(origin, dir, plane, dim, startsInside) {
        if (startsInside === void 0) { startsInside = false; }
        // assuming vectors are all normalized
        var denom = dot(plane, dir, dim);
        var abs = Math.abs(denom);
        // not parallel + must be traveling in the same direction if it starts inside
        if (abs > 0.0 && (!startsInside || denom > 0.0)) {
            return -signedDistToPlane(origin, plane, dim) / denom;
        }
        return -1;
    }
    /**
     * Calculates the adjoint of a matrix
     *
     * @ignore
     */
    function adjointMatrix(mtx, tgt, dim) {
        if (dim === 1)
            return [[1]];
        // temp is used to store cofactors of A[][]
        var sign = 1;
        var cof = getSquareMatrix(dim);
        for (var i = 0; i < dim; ++i) {
            for (var j = 0; j < dim; ++j) {
                cofactor(mtx, cof, i, j, dim);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                tgt[j][i] = (sign) * (det(cof, dim - 1));
            }
        }
        return tgt;
    }
    /**
     * Calculates the inverse of a matrix
     * @ignore
     */
    function invertMatrix(mtx, dim) {
        // custom cases for efficiency
        if (dim === 2) {
            var m00 = mtx[0][0], m01 = mtx[0][1];
            var m10 = mtx[1][0], m11 = mtx[1][1];
            var determinant = m00 * m11 - m01 * m10;
            if (determinant === 0.0)
                return null;
            var rcpDet = 1.0 / determinant;
            mtx[0][0] = m11 * rcpDet;
            mtx[0][1] = -m01 * rcpDet;
            mtx[1][0] = -m10 * rcpDet;
            mtx[1][1] = m00 * rcpDet;
        }
        else if (dim === 3) {
            var m00 = mtx[0][0], m01 = mtx[0][1], m02 = mtx[0][2];
            var m10 = mtx[1][0], m11 = mtx[1][1], m12 = mtx[1][2];
            var m20 = mtx[2][0], m21 = mtx[2][1], m22 = mtx[2][2];
            var determinant = m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
            if (determinant === 0.0)
                return null;
            var rcpDet = 1.0 / determinant;
            mtx[0][0] = (m11 * m22 - m12 * m21) * rcpDet;
            mtx[0][1] = (m02 * m21 - m01 * m22) * rcpDet;
            mtx[0][2] = (m01 * m12 - m02 * m11) * rcpDet;
            mtx[1][0] = (m12 * m20 - m10 * m22) * rcpDet;
            mtx[1][1] = (m00 * m22 - m02 * m20) * rcpDet;
            mtx[1][2] = (m02 * m10 - m00 * m12) * rcpDet;
            mtx[2][0] = (m10 * m21 - m11 * m20) * rcpDet;
            mtx[2][1] = (m01 * m20 - m00 * m21) * rcpDet;
            mtx[2][2] = (m00 * m11 - m01 * m10) * rcpDet;
            return this;
        }
        else {
            // There are faster ways of doing this, but for now, it'll do
            var determinant = det(mtx, dim);
            if (determinant === 0) {
                return null;
            }
            var rcpDet = 1.0 / determinant;
            // Find adjoint
            var adj = adjointMatrix(mtx, getSquareMatrix(dim), dim);
            // inverse(A) = adjoint/determinant
            for (var i = 0; i < dim; ++i)
                for (var j = 0; j < dim; j++)
                    mtx[i][j] = adj[i][j] * rcpDet;
        }
        return mtx;
    }
    /**
     * @ignore
     */
    function transformVector(mtx, p, tgt, dim) {
        for (var i = 0; i < dim; ++i) {
            tgt[i] = 0;
            for (var j = 0; j < dim; ++j) {
                tgt[i] += p[j] * mtx[i][j];
            }
        }
        return tgt;
    }

    /**
     * Base geometry elements.
     *
     * @author derschmale <http://www.derschmale.com>
     */
    /**
     * Ridge is a face of a facet, represented as the set of vertices forming the ridge. A line facet contains 2 ridges
     * each containing a single point. A triangle facet contains 3 ridges containing 2 vertices (segment end points), a
     * tetrahedral facet contains 4 face ridges (each containing the three triangle vertices).
     */
    var Ridge = /** @class */ (function () {
        /**
         * Creates a new ridge belonging to a facet
         */
        function Ridge(facet) {
            /**
             * The vertices of the ridge. These are always the points forming the ridge. Represented as indices into an external
             * point array.
             */
            this.verts = [];
            this.facet = facet;
        }
        /**
         * The plane containing the ridge. Once created it remains cached.
         */
        Ridge.prototype.getPlane = function (points, centroid) {
            if (!this._plane) {
                this._plane = hyperplaneFromPoints(this.verts.map(function (i) { return points[i]; }));
                if (signedDistToPlane(centroid, this._plane))
                    negate(this._plane);
            }
            return this._plane;
        };
        return Ridge;
    }());
    /**
     * In N dimensions, a facet forms an N-1 "polygon" which can be combined into an N-dimensional shape such as a
     * simplex, a convex hull, a triangulation, ...
     * A facet consists out of ridges and vertices. Care has to be taken during construction that the vertices and
     * ridges are in consistent order. The ridge's vertices must be a looping slice of (size-1) of the facet's vertices
     * starting at its corresponding index, ie:
     * - Every ridge at index I must start with the corresponding vertex at index I.
     * - If this ridge at index I stores a vertex at index N, the facet's vertex index must be (N + I) % numVerts.
     */
    var Facet = /** @class */ (function () {
        function Facet() {
            /**
             * The set of ridges for the facet. Ridges are a dimension lower than the facet (ie: points for lines, edges for
             * triangles, faces for tetrahedrons).
             */
            this.ridges = [];
            /**
             * The vertices of the facet. Represented as indices into an external point array. While they're also contained
             * in the ridges, it's useful for calculating barycentric coordinates for a facet.
             */
            this.verts = [];
        }
        return Facet;
    }());

    /**
     * Some shape construction and query code used internally.
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
    function findNeighbor(facet, ridge, facets) {
        // this simply checks if all vertices are shared for all provided facets. Not very efficient, there's probably
        // better ways to do a neighbour search. However, this is generally only applied to relatively small sets (newly
        // constructed faces).
        var src = ridge.verts;
        var numVerts = src.length;
        for (var _i = 0, facets_1 = facets; _i < facets_1.length; _i++) {
            var f = facets_1[_i];
            for (var _a = 0, _b = f.ridges; _a < _b.length; _a++) {
                var r = _b[_a];
                // do not bother if we already found them
                if (r.neighbor)
                    continue;
                var found = true;
                for (var i = 0; i < numVerts; ++i)
                    found = found && r.verts.indexOf(src[i]) >= 0;
                // all vertices are shared
                if (found) {
                    ridge.neighbor = r;
                    r.neighbor = ridge;
                    return;
                }
            }
        }
    }
    /**
     * Generates a plane for a facet. Used when constructing new facets.
     * @param facet The facet to construct a new plane for.
     * @param points The general set of points the indices refer to.
     * @param centroid An optional point that's guaranteed to be behind the plane. This can be used to orient the face
     * if the vertex order is inconsistent.
     *
     * @ignore
     */
    function generateFacetPlane(facet, points, dim, centroid) {
        var verts = facet.ridges.map(function (r) { return points[r.verts[0]]; });
        var plane = facet.plane = hyperplaneFromPoints(verts);
        if (centroid && signedDistToPlane(centroid, plane, dim) > 0.0) {
            negate(plane);
            // flip ridges for consistency
            facet.verts.reverse();
            facet.ridges.reverse();
            for (var _i = 0, _a = facet.ridges; _i < _a.length; _i++) {
                var r = _a[_i];
                r.verts.reverse();
            }
        }
    }
    /**
     * Generates the ridges for a facet based on the vertices contained in it. If ridges are already present, they're
     * considered already built and only the missing ones will be added (used only when connecting ridges). In this
     * case, the order of vertices MUST match that of the present ridges.
     *
     * @ignore
     */
    function buildRidges(facet, facets, dim) {
        var verts = facet.verts;
        var ridges = facet.ridges;
        for (var r = ridges.length; r < dim; ++r) {
            var ridge = new Ridge(facet);
            for (var v = 0; v < dim - 1; ++v) {
                ridge.verts[v] = verts[(r + v) % dim];
            }
            ridge.opposite = verts[(r + dim) % dim];
            findNeighbor(facet, ridge, facets);
            facet.ridges.push(ridge);
        }
    }
    /**
     * Combines a ridge and a point into a new facet.
     * @param ridge The ridge to extend.
     * @param p The index of the new point to build the missing ridges from.
     * @param points The array containing the point values.
     * @param facets The other facets in the shape, used to find neighbours. Usually, when constructing closed shapes,
     * these are only the newly constructed faces.
     * @param insidePoint A point guaranteed to be on the negative side of the facet plane.
     * @param dim The dimension of the facet.
     *
     * @ignore
     */
    function extendRidge(ridge, p, points, facets, insidePoint, dim) {
        // in 2D, we simply need to create 1 new facet (line) from old ridge to p
        var facet = new Facet();
        // collect all verts for this facet, which is the horizon ridge + this point
        facet.verts = ridge.verts.concat([p]);
        // the horizon ridge is part of the new facet, and gets to keep its neighbor
        facet.ridges.push(ridge);
        ridge.facet = facet;
        // the new opposite point
        ridge.opposite = p;
        buildRidges(facet, facets, dim);
        generateFacetPlane(facet, points, dim, insidePoint);
        return facet;
    }
    /**
     * Calculates the centroid ("average") for a collection of points.
     *
     * @param points The array containing all point coordinates.
     * @param indices The indices of the points for which to calculate the averages. If not provided, the first N+1
     * (simplex) points are used from the points array.
     * @param target An optional target to store the centroid in.
     * @ignore
     */
    function createCentroid(points, indices, target) {
        var index = indices ? indices[0] : 0;
        var p0 = points[index];
        var dim = p0.length;
        var numPoints = indices ? indices.length : dim + 1;
        // a point that will be internal from the very first simplex. Used to correctly orient new planes
        if (target) {
            for (var j = 0; j < dim; ++j)
                target[j] = p0[j];
        }
        else {
            target = p0.slice();
        }
        for (var j = 0; j < dim; ++j) {
            for (var i = 1; i < numPoints; ++i) {
                var index_1 = indices ? indices[i] : i;
                target[j] += points[index_1][j];
            }
            target[j] /= numPoints;
        }
        return target;
    }

    /**
     * Creates an N-simplex from N+1 points. The dimension of the points are used to
     * define the dimension of the simplex.
     *
     * @param points An array of points. Only the first N+1 points are used.
     * @param indices An optional array of indices into points to define which points in the set are used.
     *
     * @author derschmale <http://www.derschmale.com>
     */
    function createSimplex(points, indices) {
        var dim = points[0].length;
        var numVerts = dim + 1;
        var facets = [];
        for (var i = 0; i <= dim; ++i) {
            var f = new Facet();
            // collect all verts for this facet
            // the facet is made up of dim + 1 points, so cycle through these
            for (var v = 0; v < dim; ++v) {
                var index_1 = (i + v) % numVerts;
                f.verts[v] = indices ? indices[index_1] : index_1;
            }
            buildRidges(f, facets, dim);
            // an opposing face to ensure correct direction
            var index = (i + dim) % (dim + 1);
            var opposing = points[indices ? indices[index] : index];
            generateFacetPlane(f, points, dim, opposing);
            facets.push(f);
        }
        return facets;
    }

    /**
     * Randomizes the order of the elements in the array.
     */
    function shuffle(array) {
        var currentIndex = array.length, temporaryValue, randomIndex;
        // While there remain elements to shuffle...
        while (0 !== currentIndex) {
            // Pick a remaining element...
            randomIndex = Math.floor(Math.random() * currentIndex);
            currentIndex -= 1;
            // And swap it with the current element.
            temporaryValue = array[currentIndex];
            array[currentIndex] = array[randomIndex];
            array[randomIndex] = temporaryValue;
        }
        return array;
    }

    /**
     * Removes an item by removing the last item and inserting it in its place. Should only be used when the order of
     * elements is not important.
     */
    function removeIndexOutOfOrder(target, index) {
        console.assert(index < target.length, "Index out of bounds!");
        var last = target.pop();
        if (target.length > 0 && last !== target[index])
            target[index] = last;
        return index;
    }
    /**
     * Removes multiple indices out of order.
     */
    function removeIndicesOutOfOrder(target, indices) {
        // sorting them in descending order makes them easy to delete optimally. Make sure not to modify the original
        // array, we don't know if it's still used.
        var sorted = indices.slice().sort(function (a, b) { return b - a; });
        for (var _i = 0, sorted_1 = sorted; _i < sorted_1.length; _i++) {
            var i = sorted_1[_i];
            removeIndexOutOfOrder(target, i);
        }
        return target;
    }
    /**
     * Unsafely removes an item by removing the last item and inserting it in its place. Should only be used when the
     * order of elements is not important. It returns the index of the object that was removed.
     */
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

    /**
     * @ignore
     */
    var EPSILON = 0.0001;

    /**
     * Some meta-data while constructing the facets
     */
    var FacetInfo = /** @class */ (function () {
        function FacetInfo() {
            this.outsideSet = [];
            this.outsideDist = [];
        }
        return FacetInfo;
    }());
    /**
     * Gets the point in an outside set that's the furthest from the facet plane.
     *
     * @ignore
     */
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
    /**
     * Assigns all points to the outside sets of a face.
     *
     * @ignore
     */
    function generateOutsideSets(indices, points, facets, dim) {
        var outsideSets = facets.map(function (_) { return []; });
        var len = indices.length;
        for (var i = 0; i < len; ++i) {
            var index = indices[i];
            var p = points[index];
            for (var _i = 0, facets_1 = facets; _i < facets_1.length; _i++) {
                var f = facets_1[_i];
                var dist = signedDistToPlane(p, f.plane, dim);
                if (dist > EPSILON) {
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
    /**
     * Finds all faces visible to a point and their boundary ridges.
     *
     * @ignore
     */
    function getVisibleSet(p, facet, visible, horizon, dim) {
        visible.push(facet);
        facet.meta.currentPoint = p;
        for (var _i = 0, _a = facet.ridges; _i < _a.length; _i++) {
            var r = _a[_i];
            var neighbor = r.neighbor.facet;
            // already checked
            if (neighbor.meta.currentPoint === p)
                continue;
            if (signedDistToPlane(p, neighbor.plane, dim) > EPSILON)
                getVisibleSet(p, neighbor, visible, horizon, dim);
            else
                horizon.push(r);
        }
    }
    /**
     * Builds a set of new facets for a point and its horizon ridges.
     *
     * @ignore
     */
    function connectHorizonRidges(points, index, H, centroid, dim) {
        var newFacets = [];
        // link horizon ridges with new point
        for (var _i = 0, H_1 = H; _i < H_1.length; _i++) {
            var ridge = H_1[_i];
            var newFacet = extendRidge(ridge, index, points, newFacets, centroid, dim);
            newFacet.meta = new FacetInfo();
            newFacets.push(newFacet);
        }
        return newFacets;
    }
    /**
     * Returns the index with the "largest" point. The largest point is the one with the highest x coefficient, or y, z,
     * etc. if equal
     *
     * @ignore
     */
    function maximize(i, maxIndex, points, d) {
        var p = points[i];
        var max = points[maxIndex];
        for (var j = 0; j < d; ++j) {
            if (p[j] < max[j])
                return maxIndex;
            if (p[j] > max[j])
                return i;
        }
        // all the same: this only happens with duplicates, which shouldn't be in the set
        return maxIndex;
    }
    /**
     * @ignore
     */
    function minimize(i, minIndex, points, d) {
        var p = points[i];
        var min = points[minIndex];
        for (var j = 0; j < d; ++j) {
            if (p[j] > min[j])
                return minIndex;
            if (p[j] < min[j])
                return i;
        }
        // all the same: this only happens with duplicates, which shouldn't be in the set
        return minIndex;
    }
    /**
     * Tries to find the biggest shape to start with
     * @ignore
     */
    function getOptimalStart(points, d) {
        var numPoints = points.length;
        var minIndex = 0;
        var maxIndex = 0;
        // the initial axis
        for (var i = 1; i < numPoints; ++i) {
            maxIndex = maximize(i, maxIndex, points, d);
            minIndex = minimize(i, minIndex, points, d);
        }
        var indices = [minIndex, maxIndex];
        var planePts = [points[minIndex], points[maxIndex]];
        // already have 2 points, need d + 1 in total
        // in increasing dimensions, find the furthest from the current hyperplane
        for (var i = 2; i < d + 1; ++i) {
            var plane = hyperplaneFromPoints(planePts);
            var maxDist = -Infinity;
            var p = -1;
            for (var j = 0; j < numPoints; ++j) {
                var dist = Math.abs(signedDistToPlane(points[j], plane, i));
                if (dist > maxDist) {
                    maxDist = dist;
                    p = j;
                }
            }
            indices.push(p);
            planePts.push(points[p]);
        }
        return indices;
    }
    /**
     * QuickHull implements the algorithm of the same name, based on the original paper by Barber, Dobkin and Huhdanpaa.
     * We're not interested in 0- or 1-dimensional cases (the latter can simply be the extent of the point values).
     * QuickHull returns a set of indices into the original point list so we can map it to a different original ata set
     * (fe: points may be a mapping for position vectors on some scene graph object).
     *
     * @author derschmale <http://www.derschmale.com>
     */
    function quickHull(points) {
        var numPoints = points.length;
        if (numPoints === 0)
            return;
        var d = dim(points[0]);
        if (numPoints <= d) {
            throw new Error("A convex hull in " + d + " dimensions requires at least " + (d + 1) + " points.");
        }
        // initial unprocessed point indices:
        var indices = [];
        for (var i = 0; i < numPoints; ++i)
            indices.push(i);
        var simplexIndices = getOptimalStart(points, d);
        var centroid = createCentroid(points, simplexIndices);
        var facets = createSimplex(points, simplexIndices);
        for (var _i = 0, facets_2 = facets; _i < facets_2.length; _i++) {
            var f = facets_2[_i];
            f.meta = new FacetInfo();
        }
        removeIndicesOutOfOrder(indices, simplexIndices);
        shuffle(indices);
        generateOutsideSets(indices, points, facets, d);
        // do not cache facets.length, as it will keep changing
        var done = false;
        // TODO: this extra loop should not be required
        while (!done) {
            done = true;
            for (var i = 0; i < facets.length; ++i) {
                var facet = facets[i];
                var p = getFurthestPoint(facet);
                if (p !== -1) {
                    removeElementOutOfOrder(facet.meta.outsideSet, p);
                    var V = [];
                    var H = [];
                    getVisibleSet(points[p], facet, V, H, d);
                    var newFacets = connectHorizonRidges(points, p, H, centroid, d);
                    for (var _a = 0, V_1 = V; _a < V_1.length; _a++) {
                        var v = V_1[_a];
                        if (removeElementOutOfOrder(facets, v) <= i)
                            --i;
                        generateOutsideSets(v.meta.outsideSet, points, newFacets, d);
                        if (v.meta.outsideSet.length > 0)
                            done = false;
                    }
                    facets.push.apply(facets, newFacets);
                }
            }
        }
        for (var _b = 0, facets_3 = facets; _b < facets_3.length; _b++) {
            var f = facets_3[_b];
            f.meta = null;
        }
        return facets;
    }

    /**
     * 🎵🎶🎵  I lift you up to a higher dimension  🎶🎵🎶
     * 🎶🎵🎶     I'm gonna make you paraboloid     🎵🎶🎵
     *
     * @ignore
     */
    function lift(points, d) {
        var liftedDim = d + 1;
        var numPoints = points.length;
        // cache coherency
        var byteSize = liftedDim << 2; // every lifted point takes up 4 bytes per float
        var arr = new ArrayBuffer((numPoints + 1) * byteSize); // add room for one upper bound
        var i = 0;
        var bound = 0;
        var min = points[0].slice();
        var max = points[0].slice();
        var lifted = points.map(function (p) {
            var f = new Float32Array(arr, i, liftedDim);
            var s = 0;
            for (var j = 0; j < d; ++j) {
                var e = p[j];
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
        var boundPt = new Float32Array(arr, numPoints * byteSize, liftedDim);
        for (var i_1 = 0; i_1 < d; ++i_1) {
            boundPt[i_1] = (min[i_1] + max[i_1]) * .5;
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
    function delaunay(points) {
        var d = dim(points[0]);
        var numPoints = points.length;
        if (numPoints === d + 1) {
            return quickHull(points);
        }
        var lifted = lift(points, d);
        var hull = quickHull(lifted);
        return hull.filter(function (f) {
            // remove all upwards facing faces and ridges
            if (f.plane[d] >= 0.0) {
                for (var _i = 0, _a = f.ridges; _i < _a.length; _i++) {
                    var r = _a[_i];
                    if (r.neighbor)
                        r.neighbor.neighbor = null;
                }
                return false;
            }
            return true;
        });
    }

    /**
     * @ignore
     */
    var mtxCache = [];
    /**
     * @ignore
     */
    var tmpCache = [];
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
    function barycentricCoords(position, facet, points, tgt) {
        var d = dim(position);
        tgt = tgt || new Float32Array(d + 1);
        // we'll probably want to execute this function a lot of times, so let's make it efficient by not recreating these
        var mtx = mtxCache[d];
        if (!mtx) {
            mtxCache[d] = mtx = getSquareMatrix(d);
            tmpCache[d] = new Float32Array(d);
        }
        var tmp = tmpCache[d];
        var verts = facet.verts;
        var pN = points[verts[d]];
        // express all relative to pN
        for (var i = 0; i < d; ++i) {
            tmp[i] = position[i] - pN[i];
            var p = points[verts[i]];
            for (var j = 0; j < d; ++j) {
                // vectors are transposed!
                mtx[j][i] = p[j] - pN[j];
            }
        }
        invertMatrix(mtx, d);
        transformVector(mtx, tmp, tgt, d);
        tgt[d] = 1.0;
        for (var i = 0; i < d; ++i)
            tgt[d] -= tgt[i];
        return tgt;
    }

    /**
     * Provides an initial estimate to start searching, based on the facets axis-oriented bounds.
     *
     * @ignore
     */
    function findStartFacet(position, points, facets) {
        var numFacets = facets.length;
        for (var i = 0; i < numFacets; ++i) {
            var f = facets[i];
            var verts = f.verts;
            var numVerts = verts.length;
            var dim_1 = points[0].length;
            var found = true;
            for (var d = 0; d < dim_1; ++d) {
                var min = points[verts[0]][d];
                var max = min;
                // find min coordinate
                for (var v = 1; v < numVerts; ++v) {
                    var p_1 = points[verts[v]][d];
                    if (p_1 < min)
                        min = p_1;
                    else if (p_1 > max)
                        max = p_1;
                }
                var p = position[d];
                if (p < min || p >= max) {
                    found = false;
                    break;
                }
            }
            if (found)
                return f;
        }
        return null;
    }
    /**
     * Walks recursively through the neihbors of a set until the containing facet is found.
     *
     * @ignore
     */
    function walk(position, facet, points, centroid, dir, dim) {
        createCentroid(points, facet.verts, centroid);
        // so now we need to test the ray centroid -> position against the ridges and see if any intersect
        for (var i = 0; i < dim; ++i) {
            dir[i] = position[i] - centroid[i];
        }
        // using the centroid makes things easier, as the ray starting from the centroid hits the triangle face for
        // which the intersection distance is closest, so test for minT rather than doing barycentric tests.
        var hit = null;
        var minT = 1.0 - EPSILON;
        for (var _i = 0, _a = facet.ridges; _i < _a.length; _i++) {
            var r = _a[_i];
            var t = intersectRayPlane(centroid, dir, r.getPlane(points, centroid), dim, true);
            // intersection did not occur on the segment, or it's not the furthest
            if (t > EPSILON && t <= minT) {
                minT = t;
                hit = r;
            }
        }
        // if no intersection is found, we must be in the facet
        if (hit) {
            // there is an intersection, but we may have left the shape if there's no neighbour
            return hit.neighbor ?
                walk(position, hit.neighbor.facet, points, centroid, dir, dim) :
                null;
        }
        return facet;
    }
    function visibilityWalk(position, facets, points, startFacetOrEstimate) {
        var startFacet;
        if (startFacetOrEstimate === true) {
            startFacet = findStartFacet(position, points, facets);
            if (!startFacet)
                return null;
        }
        else
            startFacet = startFacetOrEstimate || facets[0];
        // this is just a single reusable object in order not to have to recreate it
        var d = dim(points[0]);
        var centroid = new Float32Array(d);
        var dir = new Float32Array(d);
        return walk(position, startFacet, points, centroid, dir, d);
    }

    exports.Facet = Facet;
    exports.Ridge = Ridge;
    exports.barycentricCoords = barycentricCoords;
    exports.createSimplex = createSimplex;
    exports.delaunay = delaunay;
    exports.quickHull = quickHull;
    exports.visibilityWalk = visibilityWalk;

    Object.defineProperty(exports, '__esModule', { value: true });

    return exports;

}({}));
