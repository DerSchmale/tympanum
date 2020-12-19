function initCanvas(canvas, drawFunc) {
    const ctx = canvas.getContext("2d");
    let w, h;
    let zOffset = 1.0;
    let rotation = 0;

    window.addEventListener("resize", function() {
        resize();
        drawFunc();
    });

    function resize()
    {
        w = canvas.width = window.innerWidth;
        h = canvas.height = window.innerHeight;
    }

    resize();

    function indexToScreen(p, points)
    {
        return toScreen(points[p]);
    }

    function toScreen(p)
    {
        const span = h * .75 * .5;
        const x = p[0];
        const y = p[1];
        let z = p.length >= 3?  p[2] : 0.0;
        let xr = x, yr = y, zr = z;

        if (rotation) {
            let c = Math.cos(rotation);
            let s = Math.sin(rotation);
            xr = c * x - s * z;
            zr = s * x + c * z;
        }

        zr += zOffset;
        // simple div Z

        return [
            w * .5 + span * xr / zr,
            h * .5 + span * yr / zr
        ]
    }

    function clear()
    {
        ctx.clearRect(0, 0, w, h);
    }

    function drawPoints(points, style, radius)
    {
        radius = radius || 2;

        ctx.fillStyle = style || "white";

        points.forEach(p => {
            ctx.beginPath();
            ctx.arc(...toScreen(p), radius, 0, Math.PI * 2);
            ctx.fill();
        })
    }

    function drawFacet(facet, points)
    {
        const numRidges = facet.ridges.length;

        for (let r = 0; r < numRidges; ++r) {
            const ridge = facet.ridges[r];
            const numVerts = ridge.verts.length;

            for (let v = 0; v < numVerts; ++v) {
                const c = indexToScreen(ridge.verts[v], points);
                if (r === 0 && v === 0)
                    ctx.moveTo(...c);
                else
                    ctx.lineTo(...c);
            }

            if (numRidges > 2)
                ctx.lineTo(...indexToScreen(ridge.verts[0], points));
        }
    }

    function drawFacetNormal(facet, points, length)
    {
        const numRidges = facet.ridges.length;
        const dim = points[0].length;
        const p = new Float32Array(dim);
        length = length || 0.05;

        let c = 0;

        for (let r = 0; r < numRidges; ++r) {
            const ridge = facet.ridges[r];
            const numVerts = ridge.verts.length;
            for (let v = 0; v < numVerts; ++v) {
                const vert = points[ridge.verts[v]];
                for (let i = 0; i < dim; ++i)
                    p[i] += vert[i];
                ++c
            }
        }

        for (let i = 0; i < dim; ++i)
            p[i] /= c;

        ctx.moveTo(...toScreen(p));

        for (let i = 0; i < dim; ++i) {
            p[i] += facet.plane[i] * length;
        }

        ctx.lineTo(...toScreen(p));
    }

    function drawFacets(facets, points, strokeStyle, lineWidth)
    {
        ctx.beginPath();
        ctx.strokeStyle = strokeStyle || "green";
        ctx.lineWidth = lineWidth || 1;

        // let i = 0;
        facets.forEach(facet => {
            // ctx.beginPath();
            // setTimeout(() => {
            //     ctx.beginPath();
                drawFacet(facet, points);
            //     ctx.stroke();
            // }, i);

            // i += 1000;
        });

        ctx.stroke();
    }

    function drawFacetNormals(facets, points, length, strokeStyle, lineWidth)
    {
        ctx.beginPath();
        ctx.strokeStyle = strokeStyle || "red";
        ctx.lineWidth = lineWidth || 1;

        facets.forEach(facet => {
            drawFacetNormal(facet, points, length);
        });

        ctx.stroke();
    }

    return {
        clear,
        drawPoints,
        drawFacets,
        drawFacetNormals,
        get rotation() { return rotation; },
        set rotation(value) { rotation = value; },
        get zOffset() { return zOffset; },
        set zOffset(value) { zOffset = value; }
    }
}

