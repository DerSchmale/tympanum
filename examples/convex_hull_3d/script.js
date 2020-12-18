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
        drawing.zOffset = 1.4;
        generate();
    }

    function generate()
    {
        points = [];

        for (let i = 0; i < 5000; ++i) {

            let azim = Math.random() * Math.PI * 2.0;
            let pol = (Math.random() - .5) * Math.PI * 2.0;
            let rad = Math.random();
            points[i] = [Math.sin(pol) * Math.cos(azim) * rad, Math.sin(pol) * Math.sin(azim) * rad, Math.cos(pol) * rad];
        }

        let time = performance.now();
        hull = TYMP.quickHull(points);

        time = performance.now() - time;
        console.log("Time to generate hull: " + time.toFixed(2) + "ms");

        draw();
    }

    function draw()
    {
        drawing.clear();
        drawing.drawPoints(points);
        drawing.drawFacets(hull, points);
        drawing.drawFacetNormals(hull, points);
    }

})();