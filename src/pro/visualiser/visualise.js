
var margin = {top: 50, right: 80, bottom: 50, left: 100}
var width = 800 - margin.left - margin.right
var height = 400 - margin.top - margin.bottom;
var maxbarheight = margin.right/2

var isdrag = false;

// X axis
var x = d3.scale.linear().range([0, width]);
var xAxis = d3.svg.axis().scale(x).orient("bottom");

// Y axis
var y = d3.scale.linear().range([height, 0]);
var yAxis = d3.svg.axis().scale(y).orient("left")
.tickFormat(function (d) {
  var prefix = d3.formatPrefix(d);
  if (Math.abs(d) >= 1000.) {
    return prefix.scale(d) + prefix.symbol;
  } else {
    return d;
  }
});

// d3 line function
var line = d3.svg.line()
.x(function(d,i) { return x(i+1) })
.y(function(d,i) { return y(d) });

// This function bins the data at x breakpoints
var histogram = function(data, domain, nbins) {
  var hist = d3.layout.histogram().range(domain).frequency(false).bins(+nbins);//x.ticks(nbins));
  return hist(data);
}

function bin_data(simulations, domain, i) {
  var data = simulations.map(function(sim) {return sim[i]});
  // bin the data
  return histogram(data, domain, d3.select("#bincount").property("value"));
}

// This function adds all the plots to the page
var runmain = function(data) {
  // d3.json("run.json", function(data) {
  var plots = Array(data.length)
  data.map(function(plotdata, fileindex) {
    var nsimulations = plotdata.data.length
    var nstages = plotdata.data[0].length

    // Create new svg
    var svg = d3.select("body").select(".plots").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + 1.5*margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var  domain=[
      d3.min(plotdata.data, function(sim){return d3.min(sim)}),
      d3.max(plotdata.data, function(sim){return d3.max(sim)})
    ];
    if (plotdata.ymin != "") {
        domain[0] = parseFloat(plotdata.ymin);
    }
    if (plotdata.ymax != "") {
        domain[1] = parseFloat(plotdata.ymax);
    }
    // set domain of axis
    y.domain(domain);
    x.domain([1, nstages]);

    line.interpolate(plotdata.interpolate);
    // add the x axis
    svg.append("g").attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")").call(xAxis);

    // add the y axis
    svg.append("g").attr("class", "y axis").call(yAxis)

    // add all the simulations
    svg.append("g").attr("class", "lines")
    var sim = svg.select("g.lines").selectAll(".sim")
    .data(plotdata.data)
    .enter().append("g")
    .attr("class", "sim");

    sim.append("path")
    // name class by simulation index
    .attr("class", function(d,i) {return "line " + i})
    // add the line
    .attr("d", function(d,i) {return line(d,i); })
    // highlight simulation on mouseover
    .on("mouseover", function(d) {
      if (isdrag == false) {
        // class name
        var classstring = d3.select(this)[0][0].className.baseVal.split(" ");
        // all corresponding simulations from all plots
        elements = document.getElementsByClassName(classstring[1]);
        for (var i = 0; i < elements.length; i++) {
          // change to light blue
          elements[i].style.stroke="#009AC7";
          // change opacity
          elements[i].style.opacity=1;
          // increase width
          elements[i].style["stroke-width"]="3px";
          // bring to front
          elements[i].parentNode.parentNode.appendChild(elements[i].parentNode);
        }
      }
    })
    // restore line on mouse out
    .on("mouseout", function(d) {
      // class name
      var classstring = d3.select(this)[0][0].className.baseVal.split(" ");
      // get all corresponding lines from all plots
      elements = document.getElementsByClassName(classstring[1]);
      for (var i = 0; i < elements.length; i++) {
        // light grey
        elements[i].style.stroke="#8D9091";
        //elements[i].style.stroke="#009AC7";
        // transparent again
        elements[i].style.opacity=0.1;
        // thin again
        elements[i].style["stroke-width"]="2px";
      }
    })
    .on("mousedown", function(yvalues, lineindex) {
      if (isdrag == false) {

        var xpos = d3.mouse(this)[0];
        var stageindex =  Math.round(xpos * (yvalues.length-1) / width)
        var color = d3.scale.pow().exponent(0.25)
        .domain([0,
          Math.max(
            Math.abs(yvalues[stageindex] - d3.min(plotdata.data, function(sim) {return sim[stageindex];})),
            Math.abs(yvalues[stageindex] - d3.max(plotdata.data, function(sim) {return sim[stageindex];}))
          )
        ])
        .range(["#8D9091", "white"]);
        var opacity = d3.scale.pow().exponent(0.25)
        .domain([0,
          Math.max(
            Math.abs(yvalues[stageindex] - d3.min(plotdata.data, function(sim) {return sim[stageindex];})),
            Math.abs(yvalues[stageindex] - d3.max(plotdata.data, function(sim) {return sim[stageindex];}))
          )
        ])
        .range([1., 0]);

        for (var i=0;i<plotdata.data.length;i++) {
          // for each line
          if (i!=lineindex) {
            // For all lines except this one, recolour
            // get all corresponding lines from all plots
            elements = document.getElementsByClassName("line " + i);
            for (var j = 0; j < elements.length; j++) {
              var diff = Math.abs(plotdata.data[i][stageindex] - yvalues[stageindex]);
              elements[j].style.stroke=color(diff);
              elements[j].style.opacity=opacity(diff);
            }
          }
        }
      }
    })
    .on("mouseup", function(yvalues, lineindex) {
      // restore colour from mousedown event
      for (var i=0;i<plotdata.data.length;i++) {
        if (i!=lineindex) {
          // get all corresponding lines from all plots
          elements = document.getElementsByClassName("line " + i);
          for (var j = 0; j < elements.length; j++) {
            // light grey
            elements[j].style.stroke="#8D9091";
            //elements[j].style.stroke="#009AC7";
            elements[j].style.opacity=0.1;
          }
        }
      }
    }); // end sim.append("path")

    // Add the text label for the x axis
    svg.append("text")
    .attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom) + ")")
    .attr("class", "axislabel")
    .text(plotdata.xlabel);

    // Add the text label for the Y axis
    svg.append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", 0 - margin.left)
    .attr("x",0 - (height / 2))
    .attr("dy", "1em")
    .attr("class", "axislabel")
    .text(plotdata.ylabel);

    svg.append("text")
    .attr("x", (width / 2))
    .attr("y", 0 - (margin.top / 2))
    .attr("text-anchor", "middle")
    .attr("class", "title")
    .text(plotdata.title);

    // add the vertical line
    svg.append("line")
    .attr("x1", width).attr("y1", 0)
    .attr("x2", width).attr("y2", height)
    .attr("class", "vertline");

    plots[fileindex] = {svg:svg, domain:domain, file:fileindex, data:plotdata.data}
  }); // end data.map

  // when the input range changes update the circle
  d3.select("#bincount").on("input", function() {
    update(+this.value);
  });

  // Initial starting radius of the circle
  update(10);

  // update the elements
  function update(bincount) {

    // adjust the text on the range slider
    d3.select("#bincount-value").text(bincount);
    d3.select("#bincount").property("value", bincount);
    if (document.getElementById("chk_histogram").checked) {
      elements = document.getElementsByClassName("vertline");
      var nstages = plots[0].data[0].length
      tocolumn = elements[0].getAttribute("x1") * (nstages-1) / width;
      plots.map(function(plot) {
        plot.svg.selectAll(".bars").remove();
      });
      addhistograms(plots);
      shifthistogram(plots, tocolumn);
    }
  }
  addhistograms(plots);
  document.getElementById("chk_histogram").addEventListener('click', function() {
    if (document.getElementById("chk_histogram").checked) {
      addhistograms(plots);
      var elements = document.getElementsByClassName("vertline");
      for (var i = 0; i < elements.length; i++) {
        elements[i].setAttribute("x1", +width);
        elements[i].setAttribute("x2", +width);
      }; // end for
    } else {
      shifthistogram(plots, 0);
      plots.map(function(plot) {
        plot.svg.selectAll("g.bars").remove();
      });
    }
  });
  // }); // end d3.json
};  // end runmain

var addhistograms = function(plots) {
  plots.map(function(plot) {
    var nstages = plot.data[0].length
    var histdata = []
    for (var i=0;i<nstages;i++) {
      histdata.push(bin_data(plot.data, plot.domain, i));
    }
    plot.histdata = histdata
    addhistogram(plots, plot)
    //plot.svg.on("dblclick", function() {
    //  // get the x location
    //  xpos = d3.mouse(this)[0];
    //  // check its within the x domain
    //  if (document.getElementById("chk_histogram").checked && xpos >= x(1) && xpos < x(nstages+0.49)) {
    //    // bin to integer weeks
    //    var tocolumn =  Math.round(xpos * (nstages-1) / width)
    //  shifthistogram(plots, tocolumn);
    // }
    //})// end plot.svg.on
  }); // end plots.map
};

var shifthistogram = function(plots, tocolumn) {
  var nstages = plots[0].histdata.length

  // bin to integer weeks
  xact = tocolumn * width / (nstages-1)
  // for all the lines
  var elements = document.getElementsByClassName("vertline");
  for (var i = 0; i < elements.length; i++) {
    elements[i].setAttribute("x1", +xact);
    elements[i].setAttribute("x2", +xact);
  }; // end for

  // rebuild histograms in correct place
  plots.map(function(plot) {
    var x = d3.scale.linear().domain(plot.domain).range([0, height]);

    // revise scale
    var y = d3.scale.linear()
    .domain([0, d3.max(plot.histdata[tocolumn], function(d){return d.y})])
    .range([maxbarheight, 0]);

    // for each bin
    for (var i = 0; i < plot.histdata[tocolumn].length; i++) {
      elements = document.getElementsByClassName("chart" + plot.file +i);
      if (elements.length > 0) {
        // update x location
        elements[0].setAttribute("x", width * (1 - tocolumn/(nstages-1)) + y(plot.histdata[tocolumn][i].y));
        // update y location
        elements[0].setAttribute("y", x(plot.histdata[tocolumn][i].x));
        // update bar width
        elements[0].setAttribute("width", margin.right/2 - y(plot.histdata[tocolumn][i].y));
      } // end if
    } // end for
  }); // end plots.map
}; // end shifthistogram

// This function initialises the histograms
function addhistogram(plots, p) {
  // x scale (corresponding to y scale of line plot)
  var x = d3.scale.linear().domain(p.domain).range([0, height]);

  // height of histogram (width in line plot)
  var y = d3.scale.linear()
  .domain([0, d3.max(p.histdata[p.histdata.length-1], function(d){return d.y})])
  .range([maxbarheight, 0]);

  // add the bars to correct svg
  var barheight = 0.
  var h;
  for (var i=0;i<p.histdata.length;i++) {
    for (var j=0;j<p.histdata[i].length;j++) {
      h = x(p.histdata[i][j].dx + p.histdata[i][j].x) - x(p.histdata[i][j].x)
      if (h > barheight) {
        barheight = h
      }
    }
  }
  p.svg.append("g").attr("class", "bars")
  var mybars = p.svg.select("g.bars").selectAll(".bar")
  .data(p.histdata[p.histdata.length-1])
  .enter();

  var xinit;
  var drag = d3.behavior.drag()
  .on("dragstart", function(d) {
    xinit = d3.mouse(this)[0]
    isdrag = true
  })
  .on("drag", function(d) {
    dx = xinit - d3.mouse(this)[0]
    xinit = d3.mouse(this)[0]
    var nstages = p.histdata.length
    var elements = document.getElementsByClassName("vertline");
    var oldx = +elements[0].getAttribute("x1");

    var tocolumn = Math.round((oldx+dx) * (nstages-1) / width)
    if (tocolumn >= nstages) {
      tocolumn = nstages - 1;
    } else if (tocolumn < 0) {
      tocolumn = 0;
    }
    shifthistogram(plots, tocolumn)
  })
  .on("dragend", function(d) {
    isdrag = false
  });

  mybars.append("rect", ".axis")
  // name class by simulation index and filename
  .attr("class", function(d,i) {return "bar chart" + p.file + i})
  // transform to correct location
  .attr("transform", function(d) { return "rotate(180) translate(" + -(width+maxbarheight) + "," + -height + ")"; })
  // y location of bar
  .attr("y", function(d) { return x(d.x);})
  // x location of bar
  .attr("x", function(d) { return y(d.y); })
  // height of the bars
  .attr("height", barheight-1)
  // width of the bars
  .attr("width", function(d) { return maxbarheight - y(d.y); })
  .call(drag);
};
