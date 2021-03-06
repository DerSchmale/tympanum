import {quickHull} from "../../build/tympanum.module.js";

(() => {
    window.onload = init;
    document.addEventListener("click", generate);

    let drawing;
    let hull;
    let points;

    function init()
    {
        drawing = initCanvas(document.getElementById("canvas2D"), draw);
        generate();
    }

    function generate()
    {
        points = [];

        /*for (let i = 0; i < 1000; ++i) {
            let ang = Math.random() * Math.PI * 2.0;
            let rad = Math.random();
            rad *= Math.sin(ang / 10.0);
            rad *= rad * 6;
            points[i] = [Math.cos(ang) * rad - 0.5, Math.sin(ang) * rad + 0.3];
        }*/

        for (let i = 0; i < 500; ++i) {

            points[i] = [
                (Math.random() - 0.5) * 2.0,
                (Math.random() - 0.5) * 2.0
            ];
        }

        let time = performance.now();
        hull = quickHull(points);

        time = performance.now() - time;
        console.log("Time to generate hull: " + time.toFixed(2) + "ms");

        draw();
    }

    function draw()
    {
        drawing.clear();
        drawing.drawFacets(hull, points);
        drawing.drawFacetNormals(hull, points);
        drawing.drawPoints(points);
    }
})();