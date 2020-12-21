import {delaunay, barycentricCoords, visibilityWalk} from "../../build/tympanum.module.js";

(() => {
    window.onload = init;
    document.addEventListener("click", reconstruct);

    let drawing;
    let triangulation;
    let points;

    function init()
    {
        drawing = initCanvas(document.getElementById("canvas2D"), draw);
        generate();
    }

    function generate()
    {
        points = [];

        for (let i = 0; i < 200; ++i) {
            points[i] = [
                (Math.random() - 0.5) * 2.0,
                (Math.random() - 0.5) * 2.0
            ];
        }

        let time = performance.now();
        triangulation = delaunay(points);

        time = performance.now() - time;
        console.log("Time to generate triangulation: " + time.toFixed(2) + "ms");

        draw();
    }

    function draw()
    {
        drawing.clear();
        drawing.drawFacets(triangulation, points);
        drawing.drawPoints(points);
    }

    function reconstruct(event)
    {
        const pos = drawing.screenTo2D(event.clientX, event.clientY);
        let time = performance.now();
        const facet = visibilityWalk(pos, triangulation, points);

        if (facet) {
            const bary = barycentricCoords(pos, facet, points);
            const p = [ 0, 0 ];

            for (let i = 0; i < bary.length; ++i) {
                const v = points[facet.verts[i]];
                const b = bary[i];
                for (let d = 0; d < 2; ++d) {
                    p[d] += v[d] * b;
                }
            }

            drawing.drawPoints([p], "blue");
        }

        time = performance.now() - time;

        console.log("Time to calculate: " + time.toFixed(2) + "ms");
    }
})();