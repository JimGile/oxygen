var xmax = 100; 
function set_x_axis_units(theseunits) {
    gasunits = theseunits; 
    if (theseunits == 'mmHg') { console.log('mmHg') }
    else if (theseunits == 'kPa') { console.log('kPa') }
    satschart.data.datasets[1].data = satline(1.0286208376646786, gasunits); 
    satschart.data.datasets[0].data = satline(); 
    satschart.options.scales.xAxes[0].scaleLabel.labelString = 'PaO2 (' + gasunits + ')'; 
    satschart.update()
}
function sat(thispo2, this_shift = 1) {
    thispo2 *= this_shift; 
    var sat = 1 / (1 + 23400 * (Math.pow((Math.pow(thispo2, 3) + 150 * thispo2), -1)))
    return sat * 100
}
function satline(s = 1, u = 'mmHg') {
    var i; el = Array(); 
    if (u == 'mmHg') {
        for (i = 0; i < xmax; i++) {
            e = sat(i, s)
            el.push({ x: i, y: e });
        }
    } else if (u == 'kPa') {
        for (i = 0; i < (10 * xmax * 101.325 / 760); i++) {
            e = sat((760 / 101.325) * (i / 10), s)
            el.push({ x: i / 10, y: e });
        }
    }
    return el
}
var ctx = document.getElementById('satschart').getContext('2d'); 
var satschart = new Chart(
        ctx, 
        { 
            type: 'scatter', 
            data: { 
                datasets: [
                    { 
                        label: 'SaO2', 
                        fill: false, 
                        borderColor: 'rgb(255, 99, 132)', 
                        showLine: true, 
                        data: satline(1.0286208376646786, gasunits), 
                    }, 
                    { 
                        label: 'normal SaO2', 
                        fill: false, 
                        borderColor: 'grey', 
                        showLine: true, 
                        data: satline(1.0286208376646786, gasunits), 
                    }
                ] 
            }, 
            options: { 
                legend: { display: true, position: 'right', }, 
                scales: { 
                    xAxes: [{ ticks: { precision: 0 }, scaleLabel: { display: true, labelString: 'PaO2 (' + gasunits + ')', fontColor: 'gray', padding: 0, } }], 
                    yAxes: [{ ticks: { suggestedMin: 0, suggestedMax: 100 }, scaleLabel: { display: true, labelString: '%HbO2', fontColor: 'gray', padding: 0 } }] 
                } 
            } 
        }
    );