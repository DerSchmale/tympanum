import {delaunay} from "../../build/tympanum.module.js";

(() => {
    window.onload = init;
    document.addEventListener("click", generate);

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
})();