//  Copyright 2017-20, Oscar Dowson.
//  This Source Code Form is subject to the terms of the Mozilla Public License,
//  v. 2.0. If a copy of the MPL was not distributed with this file, You can
//  obtain one at http://mozilla.org/MPL/2.0/.

// Compute the y domain for `plot_data`. The minimum y value is set to either
// the minimum observed value, or plot_data["ymin"] if given. A similar
// calculation is done for the maximum value ("ymax").
function compute_y_domain(plot_data) {
    domain = [
        d3.min(plot_data.data, function(sim) {
            return d3.min(sim)
        }),
        d3.max(plot_data.data, function(sim) {
            return d3.max(sim)
        })
    ];
    if (plot_data.ymin != "") {
        domain[0] = parseFloat(plot_data.ymin);
    }
    if (plot_data.ymax != "") {
        domain[1] = parseFloat(plot_data.ymax);
    }
    return domain
}

// Get the path id from an item. This is stored in the class name.
function get_path_id(item) {
    return d3.select(item)[0][0].className.baseVal.split(" ")[1]
}

// Set the stroke color, opacity, and stroke width of all lines with `path_id`.
function set_path_style(path_id, stroke, opacity, stroke_width) {
    paths = document.getElementsByClassName(path_id);
    for (i = 0; i < paths.length; i++) {
        paths[i].style.stroke = stroke;
        paths[i].style.opacity = opacity;
        paths[i].style["stroke-width"] = stroke_width;
    }
    return
}

// Bring `path_id` to the front.
function bring_to_front(path_id) {
    paths = document.getElementsByClassName(path_id);
    for (var i = 0; i < paths.length; i++) {
        paths[i].parentNode.parentNode.appendChild(paths[i].parentNode);
    }
    return
}

function main(input_data) {
    var plots = []

    var margin = {
        top: 30,
        right: 30,
        bottom: 50,
        left: 100
    }
    var width = 600 - margin.left - margin.right
    var height = 300 - margin.top - margin.bottom;

    // Cache the domain. We record `plot_domains[plot_index][stage] = {min, max}`
    var plot_domains = []
    input_data.map(function(plot_data) {
        stage_domains = []
        for (stage = 0; stage < plot_data.data[0].length; stage++) {
            stage_domains.push({
                "min": d3.min(plot_data.data, function(sim) {
                    return sim[stage]
                }),
                "max": d3.max(plot_data.data, function(sim) {
                    return sim[stage]
                })
            })
        };
        plot_domains.push(stage_domains)
    });

    input_data.map(function(plot_data, plot_index) {
        var svg = d3.select("body").select(".plots").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + 1.5 * margin.bottom)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        // Initialize the x and y scales.
        var x_scale = d3.scale.linear().range([0, width])
            .domain([1, plot_data.data[0].length]);
        var y_domain = compute_y_domain(plot_data);
        var y_scale = d3.scale.linear().range([height, 0]).domain(domain);

        // Add the x-axis.
        svg.append("g").attr("class", "x axis")
            .attr("transform", "translate(0," + height + ")").call(
                d3.svg.axis().scale(x_scale).orient("bottom")
            );
        // Add the y-axis.
        svg.append("g").attr("class", "y axis").call(
            d3.svg.axis().scale(y_scale).orient("left")
        );
        // Add a class for the lines.
        svg.append("g").attr("class", "lines")

        var line = d3.svg.line()
            .x(function(d, i) {
                return x_scale(i + 1)
            })
            .y(function(d, i) {
                return y_scale(d)
            })
            .interpolate(plot_data.interpolate);

        // add all the simulations
        var sim = svg.select("g.lines").selectAll(".sim").data(plot_data.data)
            .enter().append("g").attr("class", "sim");

        sim.append("path")
            .attr("class", function(d, i) {
                return "line " + i
            })
            .attr("d", function(d, i) {
                return line(d, i);
            })
            .on("mouseover", function(d) {
                // Highlight `this` path and bring it to the front.
                set_path_style(get_path_id(this), "#009AC7", 1, "3px");
                bring_to_front(get_path_id(this));
            })
            .on("mouseout", function(d) {
                // Set all paths back to default after mouseover finishes.
                set_path_style(get_path_id(this), "#8D9091", 0.1, "2px");
            })
            .on("mousedown", function(ydata, this_path_index) {
                // Fade away paths as they get futher away from the mouse point.
                stage = Math.round(d3.mouse(this)[0] * (ydata.length - 1) / width)
                modifying_domain = d3.scale.pow().exponent(0.4).domain(
                    [0, Math.max(
                        Math.abs(ydata[stage] - plot_domains[plot_index][stage].min),
                        Math.abs(ydata[stage] - plot_domains[plot_index][stage].max)
                    )]
                );
                for (var i = 0; i < plot_data.data.length; i++) {
                    diff = Math.abs(plot_data.data[i][stage] - ydata[stage]);
                    if (i != this_path_index) {
                        set_path_style(i,
                            modifying_domain.range(["#8D9091", "white"])(diff),
                            modifying_domain.range([1.0, 0.0])(diff),
                            "2pt"
                        )
                    }
                }
            })
            .on("mouseup", function(d, this_path_index) {
                // Restore to default all paths except the current hover.
                for (var i = 0; i < plot_data.data.length; i++) {
                    if (i != this_path_index) {
                        set_path_style(i, "#8D9091", 0.1, "2px");
                    }
                }
            });

        // Add the x-axis label.
        svg.append("text")
            .attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom) + ")")
            .attr("class", "axislabel")
            .text(plot_data.xlabel);

        // Add the y-axis label.
        svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 0 - margin.left)
            .attr("x", 0 - (height / 2))
            .attr("dy", "1em")
            .attr("class", "axislabel")
            .text(plot_data.ylabel);

        // Add the plot title.
        svg.append("text")
            .attr("x", (width / 2))
            .attr("y", 0 - (margin.top / 2))
            .attr("text-anchor", "middle")
            .attr("class", "title")
            .text(plot_data.title);

        plots.push({
            svg: svg,
            data: plot_data.data
        })
    });
};
