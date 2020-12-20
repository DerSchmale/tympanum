import {delaunay, visibilityWalk} from "../../build/tympanum.module.js";

(() => {
    window.onload = init;
    document.addEventListener("click", testContainment);

    let drawing;
    let triangulation;
    let points;
    let prevFacet = null;


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

    function testContainment(event)
    {
        const pos = drawing.screenTo2D(event.clientX, event.clientY);

        let time = performance.now();
        let facet = visibilityWalk(pos, triangulation, points, prevFacet);
        time = performance.now() - time;

        console.log("Time to find neighbour: " + time.toFixed(2) + "ms");

        prevFacet = facet;

        if (facet)
            drawing.fillFacet(facet, points);
    }
})();