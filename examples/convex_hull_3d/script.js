(() => {
    window.onload = generate;
    window.onresize = draw;


    let hull;

    let points;
    let rotation = 0;
    let moved = false;

    document.addEventListener("mousemove", e =>
    {
        if (e.buttons & 1) {
            rotation += e.movementX / 100.0;
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

    function generate()
    {
        points = [];

        for (let i = 0; i < 1000; ++i) {

            let ang = Math.random() * Math.PI * 2.0;
            let rad = Math.random();
            points[i] = [Math.cos(ang) * rad, Math.sin(ang) * rad, (Math.random() - .5) * 2.0];
        }

        let time = performance.now();
        hull = TYMP.quickHull(points);

        time = performance.now() - time;
        console.log("Time to generate hull: " + time.toFixed(2) + "ms");

        draw();
    }

    function draw()
    {
        const canvas = document.getElementById("canvas2D");
        const ctx = canvas.getContext("2d");
        const w = canvas.width = window.innerWidth;
        const h = canvas.height = window.innerHeight;
        const span = h ;

        ctx.fillStyle = "white";

        function toScreen(x, y, z)
        {
            let c = Math.cos(rotation);
            let s = Math.sin(rotation);
            let xr = c * x - s * z;
            let yr = y;
            let zr = s * x + c * z;
            zr += 4.0;
            return [
                w * .5 + span * xr / zr,
                h * .5 + span * yr / zr
            ];
        }

        function indexToScreen(i)
        {
            return toScreen(...points[i]);
        }

        points.forEach(p => {
            ctx.beginPath();
            ctx.arc(...toScreen(...p), 2, 0, Math.PI * 2);
            ctx.fill();
        })

        ctx.strokeStyle = "green";
        ctx.lineWidth = 2;

        ctx.beginPath();
        hull.forEach(facet => {
            const first = indexToScreen(facet.verts[0]);
            ctx.moveTo(...first);

            for (let i = 1; i < facet.verts.length; ++i)
                ctx.lineTo(...indexToScreen(facet.verts[i]));

            // 3 points require closing
            if (facet.verts.length > 2)
                ctx.lineTo(...first);

        })

        ctx.stroke();

        ctx.strokeStyle = "red";
        ctx.lineWidth = 1;
        ctx.beginPath();

        hull.forEach(facet => {
            const p1 = points[facet.verts[0]];
            const p2 = points[facet.verts[1]];
            const p3 = points[facet.verts[2]];
            let px = (p1[0] + p2[0] + p3[0]) / 3;
            let py = (p1[1] + p2[1] + p3[1]) / 3;
            let pz = (p1[2] + p2[2] + p3[2]) / 3;

            ctx.moveTo(...toScreen(px, py, pz));

            px += facet.plane[0] * 0.05;
            py += facet.plane[1] * 0.05;
            pz += facet.plane[2] * 0.05;

            ctx.lineTo(...toScreen(px, py, pz));
        });

        ctx.stroke();
    }

})();