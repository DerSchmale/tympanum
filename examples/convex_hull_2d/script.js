(() => {
    window.onload = generate;
    window.onresize = draw;
    document.addEventListener("click", generate);

    let hull;
    let points;

    function generate()
    {
        points = [];

        for (let i = 0; i < 200; ++i) {

            let ang = Math.random() * Math.PI * 2.0;
            let rad = Math.random();
            // rad *= Math.sin(ang / 10.0);
            // rad *= rad;
            points[i] = [Math.cos(ang) * rad, Math.sin(ang) * rad];
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
        const span = h * .75 * .5;

        ctx.fillStyle = "white";

        function toScreen(x, y)
        {
            return [
                w * .5 + span * x,
                h * .5 + span * y
            ];
        }

        function indexToScreen(i)
        {
            return toScreen(points[i][0], points[i][1]);
        }

        points.forEach(p => {
            ctx.beginPath();
            ctx.arc(...toScreen(p[0], p[1]), 2, 0, Math.PI * 2);
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
        });
        ctx.stroke();

        ctx.strokeStyle = "red";
        ctx.lineWidth = 1;
        ctx.beginPath();

        hull.forEach(facet => {
            const p1 = points[facet.verts[0]];
            const p2 = points[facet.verts[1]];
            let px = (p1[0] + p2[0]) * .5;
            let py = (p1[1] + p2[1]) * .5;

            ctx.moveTo(...toScreen(px, py));

            px += facet.plane[0] * 0.05;
            py += facet.plane[1] * 0.05;

            ctx.lineTo(...toScreen(px, py));
        });

        ctx.stroke();

    }
})();