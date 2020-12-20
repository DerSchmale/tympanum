import {quickHull} from "../../build/tympanum.module.js";

(() => {
    window.onresize = draw;
    window.onload = init;


    let hull;
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
        drawing.zOffset = 1.0;
        generate();
    }

    function generate()
    {
        points = [];

        for (let i = 0; i < 500; ++i) {
            let x = (Math.random() - .5);
            let y = (Math.random() - .5);
            let z = (x * x + y * y);
            points[i] = [ x, y, z - 0.5 ];
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