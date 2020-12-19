import {delaunay} from "../../build/tympanum.module.js";

(() => {
    window.onresize = draw;
    window.onload = init;


    let triangulation;
    let drawing;
    let points;
    let moved = false;

    document.addEventListener("mousemove", e =>
    {
        if (e.buttons & 1) {
            drawing.rotation += e.movementX / 100.0;
            draw();
            moved = true;
        }
    });

    document.addEventListener("mouseup", e => {
        if (!moved) {
            generate();
        }
        moved = false;
    });

    function init()
    {
        drawing = initCanvas(document.getElementById("canvas2D"), draw);
        drawing.zOffset = 2.0;
        generate();
    }

    function generate()
    {
        points = [];

        for (let i = 0; i < 100; ++i) {

            points[i] = [
                (Math.random() - .5) * 2.0,
                (Math.random() - .5) * 2.0,
                (Math.random() - .5) * 2.0
            ];
        }

        let time = performance.now();
        triangulation = delaunay(points);

        time = performance.now() - time;
        console.log("Time to generate hull: " + time.toFixed(2) + "ms");

        draw();
    }

    function draw()
    {
        drawing.clear();
        drawing.drawFacets(triangulation, points);
        drawing.drawPoints(points);
    }
})();