(function() {
	var draw;
	//################1.define the first layer of three plotting#####################
	//#######################define innerbox frame property##########################
	//svg width="1000" height="860">
	//	<rect class="innerBox" x="60" y="41" height="499" width="500"></rect>
	//	<rect class="innerBox" x="60" y="620" height="200" width="900"></rect>
	//	<rect class="innerBox" x="660" y="40" height="500" width="300"></rect>
	//</svg>
	
	draw = function(data) {
	var Zscale, allgenes, altpink, axislabels, bigRad, bottom, c, checkerboard, checkerboard2, chrGap, chrLowXScale, chrXScale, chrYScale, chrindex, ci, cj, cur, curXPixel, curYPixel, darkGray, darkblue, darkgreen, darkred, draw_probe, efftip, eqtltip, fasttime, g, gene, h, i, indtip, j, jitter, jitterAmount, labelcolor, left, lightGray, m, maincolor, martip, maxlod, newg, nodig, onedig, origGeneName, p, pad, peakRad, peaks, pink, probe, probesByGene, purple, pxgXscaleA, pxgXscaleX, right, slowtime, svg, tickHeight, titlecolor, top, totalChrLength, totalh, totalw, twodig, w, xloc, yloc, _i, _j, _k, _l, _len, _len1, _len10, _len2, _len3, _len4, _len5, _len6, _len7, _len8, _len9, _m, _n, _o, _p, _q, _r, _ref, _ref1, _ref10, _ref2, _ref3, _ref4, _ref5, _ref6, _ref7, _ref8, _ref9, _s;
	d3.select("p#loading").remove();
    d3.select("div#legend").style("opacity", 1);
    d3.select("div#geneinput").style("opacity", 1);
	
	totalw = 1000;
	totalh = 860;
	chrGap = 8;
	peakRad = 2;
	bigRad = 5;
	slowtime = 1000;
    fasttime = 250;
	
	pad = {
      left: 60,
      top: 40,
      right: 40,
      bottom: 40,
      inner: 10
    };
	
	left = [60,60];
	right = [560,960];
	top = [41,620,40];
	bottom = [540,820];
	h = [499,200];
	w = [500,900];
	
	//NUMBER FORMATING////////////////// 
    nodig = d3.format(".0f");
    onedig = d3.format(".1f");
    twodig = d3.format(".2f");
	//END///////////////////////////////
	
	//SET DEFAULT COLOR/////////////////
    lightGray = d3.rgb(230, 230, 230);
    darkGray = d3.rgb(200, 200, 200);
    darkblue = "darkslateblue";
    darkgreen = "darkslateblue";
    pink = "hotpink";
    altpink = "#E9CFEC";
    purple = "#8C4374";
    darkred = "crimson";
    labelcolor = "black";
    titlecolor = "blue";
    maincolor = "darkblue";
	//END///////////////////////////////
	
	svg = d3.select("div#cistrans").append("svg").attr("width", totalw).attr("height", totalh);
	for (j in left){
		svg.append("rect").attr("x", left[j]).attr("y", top[j]).attr("height", h[j]).attr("width", w[j]).attr("class", "innerBox");
	}
	//#######################end 1###################################################
	
	//###########################2.define chessboard layout##########################
	
	checkerboard = svg.append("g").attr("id", "checkerboard");// define upper left plot as checkerboard
	checkerboard2 = svg.append("g").attr("id", "checkerboard2"); // define lower plot as checkerboard2
	
	//##########  Calculate total chromosome length (in bp) SUM	#########
    totalChrLength = 0;
    _ref = data.chrnames; //{"chr" : {"1": {"start_Mbp": 3631,"end_Mbp": 30425192,"length_bp":30427671},...}}
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      c = _ref[_i];
	  //chromosome length in cM
      data.chr[c].length_bp = data.chr[c].end_Mbp - data.chr[c].start_Mbp;
      totalChrLength += data.chr[c].length_bp;
    }
	//#################################end###############################
	
	//################################Upper left plot chess board box property############################
	
	chrXScale = {};
    chrYScale = {};
    curXPixel = left[0] + peakRad; // curXPixel = 62;   left[0]=60, peakRad = 2;
    curYPixel = bottom[0] - peakRad; // curYPixel = 542;   bottom[0]=540 peakRad = 2;
    _ref1 = data.chrnames;
    for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
      c = _ref1[_j]; //c = chrname;
	  //define the width of chr will be presented in the whole chrmosome (x-axis).
      data.chr[c].length_pixel = Math.round((w[0] - peakRad * 2) * data.chr[c].length_bp / totalChrLength); // 	data.chr[c].length_pixel = (500-2*2)* data.chr.total_length_bp/totalChrLength
      data.chr[c].start_Xpixel = curXPixel; //62
      data.chr[c].end_Xpixel = curXPixel + data.chr[c].length_pixel - 1; //61+ chr_width
      data.chr[c].start_Ypixel = curYPixel; //542
      data.chr[c].end_Ypixel = curYPixel - (data.chr[c].length_pixel - 1); //543-chr_width
	  //map from an input domain to an output range.
      chrXScale[c] = d3.scale.linear().domain([data.chr[c].start_Mbp, data.chr[c].end_Mbp]).range([data.chr[c].start_Xpixel, data.chr[c].end_Xpixel]).clamp(true); //clamp(): enable or disable clamping of the output range
      chrYScale[c] = d3.scale.linear().domain([data.chr[c].start_Mbp, data.chr[c].end_Mbp]).range([data.chr[c].start_Ypixel, data.chr[c].end_Ypixel]).clamp(true);
      curXPixel += data.chr[c].length_pixel; // move the cursor to the next chr from left to right.
      curYPixel -= data.chr[c].length_pixel; // move to cursor to the next chr from up to down
    }
	//############compensatation#################
	data.chr["1"].start_Xpixel = left[0]; //60
    data.chr["1"].start_Ypixel = bottom[0]; //540
    data.chr["5"].end_Xpixel = right[0]; //560
    data.chr["5"].end_Ypixel = top[0]; // 41
    //#####################################end#######################################


	
	//###################draw chess board layout for upper left plot#################
	_ref2 = data.chrnames;
    for (i = _n = 0, _len2 = _ref2.length; _n < _len2; i = ++_n) {
      ci = _ref2[i];
      _ref3 = data.chrnames; //chrnames_list
      for (j = _o = 0, _len3 = _ref3.length; _o < _len3; j = ++_o) {
        cj = _ref3[j];
        if ((i + j) % 2 === 0) { // if i+j is even number
		  //draw chess board with dark gray color
          checkerboard.append("rect").attr("x", data.chr[ci].start_Xpixel).attr("width", data.chr[ci].end_Xpixel - data.chr[ci].start_Xpixel).attr("y", data.chr[cj].end_Ypixel).attr("height", data.chr[cj].start_Ypixel - data.chr[cj].end_Ypixel).attr("stroke", "none").attr("fill", darkGray).style("pointer-events", "none");
        }
      }
    }
	//###################################end#########################################
   
   	//################################lower plot chess board box property############################
	chrLowXScale = {};  
	cur = Math.round(pad.left + chrGap / 2); // (60+8/2)
	_ref4 = data.chrnames;
    for (_k = 0, _len4 = _ref4.length; _k < _len4; _k++) {
      c = _ref4[_k];
      data.chr[c].start_lowerXpixel = cur; // 64
	  // 64 + round(900-8*nr_chr/total_chr_length*chr_cM)
      data.chr[c].end_lowerXpixel = cur + Math.round((w[1] - chrGap * data.chrnames.length) / totalChrLength * data.chr[c].length_bp);
      chrLowXScale[c] = d3.scale.linear().domain([data.chr[c].start_Mbp, data.chr[c].end_Mbp]).range([data.chr[c].start_lowerXpixel, data.chr[c].end_lowerXpixel]);
      cur = data.chr[c].end_lowerXpixel + chrGap; // each loop cursor move to the end position +2
    }
    //#####################################end#######################################
	
	//###################draw chess board layout for lower plot#################
	_ref5 = data.chrnames;
    for (i = _p = 0, _len5 = _ref5.length; _p < _len5; i = ++_p) {
      ci = _ref5[i];
      if (i % 2 === 0) {
        checkerboard2.append("rect").attr("x", data.chr[ci].start_lowerXpixel - chrGap / 2).attr("width", data.chr[ci].end_lowerXpixel - data.chr[ci].start_lowerXpixel + chrGap).attr("y", top[1]).attr("height", h[1]).attr("stroke", "none").attr("fill", darkGray).style("pointer-events", "none");
      }
    }
	//###################################end#########################################
	
	//#########################draw axis of upper left and lower plot###########################
	
	//################################### X axis #######################################

    axislabels = svg.append("g").attr("id", "axislabels").style("pointer-events", "none");
	//################# x axis text of upper left plot################
    axislabels.append("g").attr("id", "topleftX").selectAll("empty").data(data.chrnames).enter().append("text").text(function(d) {
      return d;
    }).attr("x", function(d) {
      return (data.chr[d].start_Xpixel + data.chr[d].end_Xpixel) / 2;
    }).attr("y", bottom[0] + pad.bottom * 0.3).attr("fill", labelcolor);
	
	//################################### END #######################################
	
	//################################### Y axis #######################################
	
	//################# y axis text of upper left plot################
    axislabels.append("g").attr("id", "topleftY").selectAll("empty").data(data.chrnames).enter().append("text").text(function(d) {
      return d;
    }).attr("x", left[0] - pad.left * 0.15).attr("y", function(d) {
      return (data.chr[d].start_Ypixel + data.chr[d].end_Ypixel) / 2;
    }).style("text-anchor", "end").attr("fill", labelcolor);
    //################# x axis text of lower plot################
	axislabels.append("g").attr("id", "bottomX").selectAll("empty").data(data.chrnames).enter().append("text").text(function(d) {
      return d;
    }).attr("x", function(d) {
      return (data.chr[d].start_lowerXpixel + data.chr[d].end_lowerXpixel) / 2;
    }).attr("y", bottom[1] + pad.bottom * 0.3).attr("fill", labelcolor);
	
	//################################### END #######################################
	
	//################################### X axis label text #######################################
	//############################### add x-axis name upper left ##################################
	axislabels.append("text").text("Marker position (cM)").attr("x", (left[0] + right[0]) / 2).attr("y", bottom[0] + pad.bottom * 0.75).attr("fill", titlecolor).attr("text-anchor", "middle");
	//############################### add x-axis name lower #######################################
    axislabels.append("text").text("Marker position (cM)").attr("x", (left[1] + right[1]) / 2).attr("y", bottom[1] + pad.bottom * 0.75).attr("fill", titlecolor).attr("text-anchor", "middle");
	
	//################################### Y axis label text #######################################
	//############################### add Y-axis name upper left ##################################
    xloc = left[0] - pad.left * 0.65; // xloc and yloc together locate the title lable in the middle of y axis.
    yloc = (top[0] + bottom[0]) / 2;
    axislabels.append("text").text("Gene position (cM)").attr("x", xloc).attr("y", yloc).attr("transform", "rotate(270," + xloc + "," + yloc + ")").style("text-anchor", "middle").attr("fill", titlecolor);
	//############################### add Y-axis name lower #######################################
    xloc = left[1] - pad.left * 0.65;
    yloc = (top[1] + bottom[1]) / 2;
    axislabels.append("text").text("LOD score").attr("x", xloc).attr("y", yloc).attr("transform", "rotate(270," + xloc + "," + yloc + ")").style("text-anchor", "middle").attr("fill", titlecolor);
	//################################### END #######################################
	
	//sort the peaks by peaks.lod in the ascending order
	//selection.sort(comparator)
	//The comparator function is passed two data elements a and b to compare. If negative, then a should be before b; if positive, then a should be after b; otherwise, a and b are considered equal and the order is arbitrary. 
    data.peaks = data.peaks.sort(function(a, b) {
      if (a.lod < b.lod) {
        return -1;
      } else {
        return +1;
      }
    });
	//sort the data.exp by peaks.lod in the ascending order
	data.exp = data.exp.sort(function(a, b) {
      if (a.lod < b.lod) {
        return -1;
      } else {
        return +1;
      }
    });
	
	//###??????????????????????????????????????? ?????????###
    Zscale = d3.scale.linear().domain([0, 25]).range([0, 1]);
	
	//Marker-Probe (eQTL) tooltip/////////
    eqtltip = d3.svg.tip().orient("right").padding(3).text(function(z) {
      return "" + z.gene + " (LOD = " + (onedig(z.lod)) + ")";
    }).attr("class", "d3-tip").attr("id", "eqtltip");
	//end/////////////////////////////////
	
	//marker tooltip//////////////////////
    martip = d3.svg.tip().orient("right").padding(3).text(function(z) {
      return z;
    }).attr("class", "d3-tip").attr("id", "martip");
	//end/////////////////////////////////
	
	
	//######################draw dot plot#########################
	draw_plot = function(selected_gene) {
	  
      var chr, curves, draw_pxgXaxis, ensembl, lastMarker, lod, lodcurve, lodcurve_yScale, markerClick, minlod,minlod_marker,maxlod_marker, meanmarks, mgi, plotPXG, pos, probeaxes, pxgXaxis, pxgYaxis, pxgYscale, revPXG, ticks, title, titletext, xlink, yaxis, _len10, _len11, _len12, _ref10, _ref11, _ref12,_ref13, _s, _t, _u;
      svg.selectAll(".probe_data").remove();
      d3.select("text#pxgtitle").text("");
      svg.selectAll(".plotPXG").remove();
	
      _ref10 = data.exp;
	  maxlod = _ref10[_ref10.length-1].lod;
	  //maxlod_marker = _ref10[_ref10.length-1].marker;
	  minlod = _ref10[0].lod;
	  //minlod_marker = _ref10[0].marker;
	  
	  //############################map y axis range to the width of the lower plot#############################
	  
      lodcurve_yScale = d3.scale.linear().domain([minlod * 1.05, maxlod * 1.05]).range([bottom[1], top[1]]);//540,620
	  //add a container to group objects, in this case add y axis of lower plot
      yaxis = svg.append("g").attr("class", "probe_data").attr("id", "loweryaxis");
	  ticks = lodcurve_yScale.ticks(6);
      yaxis.selectAll("empty").data(ticks).enter().append("line").attr("y1", function(d) {
        return lodcurve_yScale(d);
      }).attr("y2", function(d) {
        return lodcurve_yScale(d);
      }).attr("x1", left[1]).attr("x2", right[1]).attr("stroke", "white").attr("stroke-width", "1");
      yaxis.selectAll("empty").data(ticks).enter().append("text").text(function(d) {
        if (maxlod > 10) { // for the sacle of y axis of lower plot
          return nodig(d);
        } else {
          return onedig(d); // if maxlod<10: rescale y axis with one digit numbers
        }
      }).attr("y", function(d) {
        return lodcurve_yScale(d);
      }).attr("x", left[1] - pad.left * 0.1).style("text-anchor", "end").attr("fill", labelcolor);
      yaxis.append("line").attr("y1", lodcurve_yScale(5)).attr("y2", lodcurve_yScale(5)).attr("x1", left[1]).attr("x2", right[1]).attr("stroke", purple).attr("stroke-width", "1").attr("stroke-dasharray", "2,2");
      lodcurve = function(c) {
        return d3.svg.line().x(function(p) {
          return chrLowXScale[c](data.pmark[p].pos_Mbp);
        }).y(function(p) {
		  _ref6 = data.exp;
		  var p_lod = 0;
		  for(i = 0; i<_ref6.length; i++){
			if (data.exp[i].gene == selected_gene){
				if (data.exp[i].marker == p){
					p_lod = data.exp[i].lod;
				}
			}
		  }	
		  return lodcurve_yScale(p_lod); 							
        });
      };
	  //#######################################end#################################################
	  
	  
	  
	  curves = svg.append("g").attr("id", "curves").attr("class", "probe_data");
      _ref11 = data.chrnames;
      for (_t = 0, _len11 = _ref11.length; _t < _len11; _t++) {
        c = _ref11[_t];
        curves.append("path").datum(data.pmarknames[c]).attr("d", lodcurve(c)).attr("class", "thickline").attr("stroke", darkblue).style("pointer-events", "none").attr("fill", "none");
      }
	  
      
	  titletext = selected_gene;
      probeaxes = svg.append("g").attr("id", "probe_data_axes").attr("class", "probe_data");
      //gene = data.probes[probe_data.probe].gene;
      ensembl = "http://plants.ensembl.org/Arabidopsis_thaliana/Gene/Summary?g=" + selected_gene;
      mgi = "https://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&show_obsolete=F&name="+selected_gene+"&sub_type=gene&SEARCH_EXACT=4&SEARCH_CONTAINS=1";
      if (selected_gene !== null) {
        titletext += " (" + selected_gene + ")";
        xlink = probeaxes.append("a").attr("xlink:href", ensembl);
        xlink.append("text").text(titletext).attr("x", (left[1] + right[1]) / 2).attr("y", top[1] - pad.top / 2).attr("fill", maincolor).style("font-size", "18px");
      } else {
        probeaxes.append("text").text(titletext).attr("x", (left[1] + right[1]) / 2).attr("y", top[1] - pad.top / 2).attr("fill", maincolor).style("font-size", "18px");
      }
      svg.append("rect").attr("class", "probe_data").attr("x", left[1]).attr("y", top[1]).attr("height", h[1]).attr("width", w[1]).attr("class", "outerBox");
      svg.append("circle").attr("class", "probe_data").attr("id", "probe_circle").attr("cx", chrLowXScale[data.gene[selected_gene].chr](data.gene[selected_gene].pos_Mbp)).attr("cy", top[1]).attr("r", bigRad).attr("fill", pink).attr("stroke", darkblue).attr("stroke-width", 1).attr("opacity", 1);
  
	  //############# Calcuate the number of markers and save to markerClick with default value of 0 ###################
	  markerClick = {};
      _ref12 = data.markers;
      for (_u = 0, _len12 = _ref12.length; _u < _len12; _u++) {
        m = _ref12[_u];
        markerClick[m] = 0;
      }
	  //########################## end ############################
	  
      lastMarker = "";
      svg.append("g").attr("id", "markerCircle").attr("class", "probe_data").selectAll("empty").data(data.markers).enter().append("circle").attr("class", "probe_data").attr("id", function(td) {
        return "marker_" + td;
      }).attr("cx", function(td) {
        return chrLowXScale[data.pmark[td].chr](data.pmark[td].pos_Mbp);
      }).attr("cy", function(td) {
		
		_ref7 = data.exp;
		var td_lod = 0;
		for(i = 0; i<_ref7.length; i++){
		  if (data.exp[i].gene == selected_gene){
			if (data.exp[i].marker == td){
			  td_lod = data.exp[i].lod;
			}
		  }
		}	
        return lodcurve_yScale(td_lod);
      }).attr("r", bigRad).attr("fill", purple).attr("stroke", "none").attr("stroke-width", "2").attr("opacity", 0).on("mouseover", function(td) {
        if (!markerClick[td]) {
          d3.select(this).attr("opacity", 1);
        }
        return martip.call(this, td);
      }).on("mouseout", function(td) {
        d3.select(this).attr("opacity", markerClick[td]);
        return d3.selectAll("#martip").remove();
      }).on("click", function(td) {
        var chr, pos, title;
        pos = data.pmark[td].pos_cM;
        chr = data.pmark[td].chr;
        title = "" + td + " (chr " + chr + ", " + (onedig(pos)) + " cM)";
		});
		//end draw lower plot
	  };
	  //############################## upper left dot eQTL plot ##############################
	  chrindex = {};
      _ref10 = data.chrnames;
      for (i = _s = 0, _len10 = _ref10.length; _s < _len10; i = ++_s) {
        c = _ref10[i];
        chrindex[c] = i;
      }
	
      peaks = svg.append("g").attr("id", "peaks").selectAll("empty").data(data.peaks).enter().append("circle").attr("class", function(d) {
        return "probe_" + d.gene; //d--->data(data.peaks) iteration of all object in peaks.
      }).attr("cx", function(d) {
	    return chrXScale[data.pmark[d.marker].chr](data.pmark[d.marker].pos_Mbp);
      }).attr("cy", function(d) {
        return chrYScale[data.gene[d.gene].chr](data.gene[d.gene].pos_Mbp);
      }).attr("r", peakRad).attr("stroke", "none").attr("fill", function(d) {
        if (chrindex[data.pmark[d.marker].chr] % 2 === 0) {
          return darkblue;
        } else {
          return darkgreen;
        }
      }).attr("opacity", function(d) {
        return Zscale(d.lod);
      }).on("mouseover", function(d) {
        d3.selectAll("circle.probe_" + d.gene).attr("r", bigRad).attr("fill", pink).attr("stroke", darkblue).attr("stroke-width", 1).attr("opacity", 1);
        return eqtltip.call(this, d);
      }).on("mouseout", function(d) {
        d3.selectAll("circle.probe_" + d.gene).attr("r", peakRad).attr("fill", function(d) {
          if (chrindex[data.pmark[d.marker].chr] % 2 === 0) {
            return darkblue;
          } else {
            return darkgreen;
          }
        }).attr("stroke", "none").attr("opacity", function(d) {
          return Zscale(d.lod);
        });
        return d3.selectAll("#eqtltip").remove();
      });
	  //##################### eQTL onclick event#################
	  /*
	  .on("click", function(d) {
        draw_plot(d.gene);     
        $('input#genesymbol').each(function() {
          return $(this).val($(this).data('default')).addClass('inactive');
        });
        return d3.select("a#currentgenesymbol").html("");
      });*/ 
	  var _gene = "AT3G50500";
	  //var _gene = "AT5G53940";
	  //var peak_marker = "K9I9"; 
	  draw_plot(_gene);
    };
	d3.json("data/test.json",draw);
}).call(this);
  
  