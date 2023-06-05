var serverscript = "https://baillielab.roslin.ed.ac.uk/cgi-bin/oxygen_delivery.py";
var es_script = "https://baillielab.roslin.ed.ac.uk/cgi-bin/abg/shunt.py";

var gasunits = 'kPa';
var showerrors = "yes"; // send console message when submission fails

function getUrl(url) {
  // get query string from url (optional) or window
  var queryString = url ? url.split('?')[1] : window.location.search.slice(1);
  // we'll store the parameters here
  var obj = {};
  // if query string exists
  if (queryString) {
    // stuff after # is not part of query string, so get rid of it
    queryString = queryString.split('#')[0];
    // split our query string into its component parts
    var arr = queryString.split('&');
    for (var i = 0; i < arr.length; i++) {
      // separate the keys and the values
      var a = arr[i].split('=');
      // set parameter name and value (use 'true' if empty)
      var paramName = a[0];
      var paramValue = typeof (a[1]) === 'undefined' ? true : a[1];
      // (optional) keep case consistent
      paramName = paramName.toLowerCase();
      if (typeof paramValue === 'string') paramValue = paramValue.toLowerCase();
      // if the paramName ends with square brackets, e.g. colors[] or colors[2]
      if (paramName.match(/\[(\d+)?\]$/)) {
        // create key if it doesn't exist
        var key = paramName.replace(/\[(\d+)?\]/, '');
        if (!obj[key]) obj[key] = [];
        // if it's an indexed array e.g. colors[2]
        if (paramName.match(/\[\d+\]$/)) {
          // get the index value and add the entry at the appropriate position
          var index = /\[(\d+)\]/.exec(paramName)[1];
          obj[key][index] = paramValue;
        } else {
          // otherwise add the value to the end of the array
          obj[key].push(paramValue);
        }
      } else {
        // we're dealing with a string
        if (!obj[paramName]) {
          // if it doesn't exist, create property
          obj[paramName] = paramValue;
        } else if (obj[paramName] && typeof obj[paramName] === 'string'){
          // if property does exist and it's a string, convert it to an array
          obj[paramName] = [obj[paramName]];
          obj[paramName].push(paramValue);
        } else {
          // otherwise add the property
          obj[paramName].push(paramValue);
        }
      }
    }
  }
  return obj;
}

/*
from https://www.sitepoint.com/get-url-parameters-with-javascript/
getUrl().product; // 'shirt'
getUrl().color; // 'blue'
getUrl().newuser; // true
getUrl().nonexistent; // undefined
getUrl('http://test.com/?a=abc').a; // 'abc'
*/

//set starting value for global variables
var fio2 = getUrl().f || 21;
var RR = getUrl().rr || 12;
var VT = getUrl().vt || 475;
var VD =  getUrl().vd || 110;
var RQ =  getUrl().rq || 0.8;
var alt =  getUrl().a || 0;
var Temp = getUrl().t || 36.5;
var Hb = getUrl().hb || 15;
var CO =  getUrl().co || 6.5;
var BE =  getUrl().be || 0;
var DPG = getUrl().dpg || 4.65;
var pulm_shunt = getUrl().sh || 0.03;
var tissue_shunt = getUrl().tsh || 0.05;
var Vc = getUrl().vc || 75;
var DmO2 = getUrl().dm || 300;
var sigma = getUrl().sigma || 0;
var VO2 = getUrl().vo || 250;
//view
var currentgraph = getUrl().currentgraph || 'c'; // "c", "s", "t", "v", "a"
var advm = getUrl().advm || 'n';  //"n" - normal, "a" - advanced
var nums = getUrl().nums || 0; // 0 or 1

document.ajax.fio2.value = fio2;
document.ajax.RR.value = RR;
document.ajax.VT.value = VT;
document.ajax.VD.value = VD;
document.ajax.RQ.value = RQ;
document.ajax.alt.value = alt;
document.ajax.Temp.value = Temp;
document.ajax.Hb.value = Hb;
document.ajax.CO.value = CO;
document.ajax.BE.value = BE;
document.ajax.DPG.value = DPG;
document.ajax.pulm_shunt.value = pulm_shunt;
document.ajax.tissue_shunt.value = tissue_shunt;
document.ajax.Vc.value = Vc;
document.ajax.DmO2.value = DmO2;
document.ajax.VO2.value = VO2;

var inputboxes={'RR':RR, 'VT':VT, 'VD':VD, 'fio2':fio2, 'alt':alt, 'CO':CO, 'pulm_shunt':pulm_shunt, 'Vc':Vc, 'DmO2':DmO2, 'sigma':sigma, 'Hb':Hb, 'BE':BE, 'DPG':DPG, 'Hct':0.35, 'VO2':VO2, 'Temp':Temp, 'tissue_shunt':tissue_shunt, 'RQ':RQ};
var returned_numbers=[21.26 , 19.96 , 13.6823308318 , 13.6823450548 , 13.0792809165 , 6.27708528477 , 6.38081851231 , 13.1 , 5.3 , 7.4 , 39.77 , 24.4 , 97.3];

function createInstance(){
        var req = null;
        if (window.XMLHttpRequest)
        {req = new XMLHttpRequest();}
        else if (window.ActiveXObject) {
        try { req = new ActiveXObject("Msxml2.XMLHTTP");}
        catch (e) { try { req = new ActiveXObject("Microsoft.XMLHTTP");}
        catch (e) { alert("XHR not created");}}}
        return req;
    };

function working(data){
        var element = document.getElementById('working');
        element.innerHTML = data;
    }

function submitForm(mode){
    if (mode=="initiate"){working("Loading model...")}
    else {working("Calculating...")}
    //shade displayed values grey
    // enter code here...

    var req =  createInstance();
    var advm = document.getElementById('versionbox').checked;

    fio2 = Number(document.ajax.fio2.value);
    RR = Number(document.ajax.RR.value);
    VT = Number(document.ajax.VT.value);
    VD = Number(document.ajax.VD.value);
    alt = Number(document.ajax.alt.value);
    CO = Number(document.ajax.CO.value);
    pulm_shunt = Number(document.ajax.pulm_shunt.value);
    DmO2 = Number(document.ajax.DmO2.value);
    Vc = Number(document.ajax.Vc.value);
    Hb = Number(document.ajax.Hb.value);
    BE = Number(document.ajax.BE.value);
    DPG = Number(document.ajax.DPG.value);
    VO2 = Number(document.ajax.VO2.value);
    Temp = Number(document.ajax.Temp.value);
    tissue_shunt = Number(document.ajax.tissue_shunt.value);
    RQ = Number(document.ajax.RQ.value);
    Hct = Number(document.ajax.Hct.value);
    //sigma = Number(document.ajax.sigma.value);
    VT_unit = document.ajax.VT_unit.value;
    VD_unit = document.ajax.VD_unit.value;
    alt_unit = document.ajax.alt_unit.value;
    Temp_unit = document.ajax.Temp_unit.value;
    Hb_unit = document.ajax.Hb_unit.value;
    DPG_unit = document.ajax.DPG_unit.value;
    CO_unit = document.ajax.CO_unit.value;
    Vc_unit = document.ajax.Vc_unit.value;

    Qecmo = Number(document.ajax.qecmo.value);
    vv_ecmo = ''
    if (document.getElementById('venacava').checked){vv_ecmo = 'venacava'};
    va_ecmo = ''
    if (document.getElementById('aorta').checked){va_ecmo = 'aorta'};
    ecmosites = vv_ecmo + '|' + va_ecmo;

    directlink = document.URL.split('?')[0] + "?f="+fio2+"&rr="+RR+"&vt="+VT+"&vd="+VD+"&a="+alt+"&co="+CO+"&sh="+pulm_shunt+"&dm="+DmO2+"&vc="+Vc+"&hb="+Hb+"&be="+BE+"&dpg="+DPG+"&vo="+VO2+"&t="+Temp+"&tsh="+tissue_shunt+"&rq="+RQ+"&currentgraph="+currentgraph+"&advm="+advm+"&nums="+nums+"&sigma"+sigma

    var data =  "&advm="+advm+
                "&fio2="+fio2+
                "&RR="+RR+
                "&VT="+VT+
                "&VT_unit="+VT_unit+
                "&VD="+VD+
                "&VD_unit="+VD_unit+
                "&RQ="+RQ+
                "&alt="+alt+
                "&alt_unit="+alt_unit+
                "&Temp="+Temp+
                "&Temp_unit="+Temp_unit+
                "&Hb="+Hb+
                "&Hb_unit="+Hb_unit+
                "&Hct="+Hct+
                "&BE="+BE+
                "&DPG="+DPG+
                "&DPG_unit="+DPG_unit+
                "&pulm_shunt="+pulm_shunt+
                "&tissue_shunt="+tissue_shunt+
                "&Vc="+Vc+
                "&Vc_unit="+Vc_unit+
                "&DmO2="+DmO2+
                "&sigma="+sigma+
                "&VO2="+VO2+
                "&CO="+CO+
                "&CO_unit="+CO_unit+
                "&Qecmo="+Qecmo+
                "&ecmosites="+ecmosites
                ;

    req.onreadystatechange = function(){
        if(req.readyState == 4){
            if(req.status == 200){
                output=req.responseText;
                console.log("output:", output)
                output = output.split(':::')[1]
                returned_numbers = output.split("|");

                results = {}
                for (index in returned_numbers){
                    s = returned_numbers[index].split(":")
                    results[s[0].replace(/\s/g, '')] = parseFloat(s[1])
                }
                update_chart(results);
                populate_abg(results);
                submitES(results);

                //remove shading from displayed values - revert to "inherit" or equivalent
                working("Done");
                }
            else
                {
                if (showerrors == "yes")
                    {
                    console.log("Error at script submission: returned status code " + req.status + " " + req.statusText +" "+ data);
                    }
                }
            }
        }
    req.open("POST", serverscript, true);
    req.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
    req.send(data);
    }

function populate_abg(res){
    document.getElementById("abg_pao2").innerHTML = res['PaO2'].toFixed(1);
    document.getElementById("abg_paco2").innerHTML = res['PaCO2'].toFixed(1);
    document.getElementById("abg_pH").innerHTML = res['pH_a'].toFixed(2);
    document.getElementById("abg_H").innerHTML = (Math.pow(10,-res['pH_a'])*1000000000).toFixed(1);
    document.getElementById("abg_bicarb").innerHTML = (res['HCO3_a']*1000).toFixed(1);
    document.getElementById("abg_sao2").innerHTML = (res['SaO2']*100).toFixed(1);
    }

function ABG(extension){
    populate_abg();
}

function advanced_toggle_1(id){
        var element = document.getElementsByTagName(id);
        for(var x = 0; x < element.length; x++){
            name = element[x].getAttribute('name');
            if (name == 'advanced'){
                if (document.getElementById('versionbox').checked == true){
                        element[x].style.display = '';
                }
                else{
                        element[x].style.display = 'none';
                }
            }
        }
    }

function advanced_toggle()
    {
    if (document.getElementById('versionbox').checked == true){submitForm("initiate")};
    advanced_toggle_1('div');
    advanced_toggle_1('dl');
    }

function getlink()
    {
        submitForm();
        window.location = directlink;
    }

function toggle(thisdiv)
    {
    if (document.getElementById(thisdiv).style.display == '')
        {document.getElementById(thisdiv).style.display = 'none';}
    else {document.getElementById(thisdiv).style.display = '';}
    }

function changed(e, name)
    {
    var element = document.getElementById(name)
    if (element.value != inputboxes[name])
        {
        element.setAttribute('class','text_input_altered');
        }
    var keycode;
    if (window.event) keycode = window.event.keyCode;
    else if (e) keycode = e.which;
    else return true;

    if (keycode == 13) //ie if the button pressed was ENTER
       {
       submitForm();
       return false;
       }
    else
       return true;
    }

function resetall()
    {
    for (var key in inputboxes){
        document.getElementById(key).setAttribute('class','text_input');
        }
    submitForm();
    }

function startup()
    {
    document.ajax.reset();
    initMove();
    submitForm('initiate');
    advanced_toggle();
    }

function abgsubmit()
    {
    var fio2 = Number(document.ajax.fio2.value);
    var abg_pao2 = Number(document.getElementById('abg_pao2').innerHTML)
    var abg_paco2 = Number(document.getElementById('abg_paco2').innerHTML)
    var abg_H = Number(document.getElementById('abg_H').innerHTML)
    window.location = '/abg/arterial_blood_gas_calculator.php?f='+fio2+'&o='+abg_pao2+'&c='+abg_paco2+'&h='+abg_H; //+"?backurl='"+directlink+"'";
    }

function submitES(res){
    document.getElementById('es').style.color = 'grey';
    var data =  "&online=yes"+
                "&fio2="+document.ajax.fio2.value/100+
                "&pao2="+res['PaO2']+
                "&paco2="+res['PaCO2']+
                "&pH="+res['pH_a']+
                "&Temp="+Temp+
                "&gasunit="+'kPa'+
                "&acidunit="+'pH'
                ;
    console.log(data)
    var req =  createInstance();
    req.onreadystatechange = function(){
            if(req.readyState == 4){
                if(req.status == 200){
                    output=req.responseText;
                    console.log("ES:", output)
                    shunt = (parseFloat(output)*100)
                    console.log("ES:", shunt)
                    shunt = shunt.toFixed(1)
                    console.log("ES:", shunt)
                    document.getElementById('es').innerHTML = shunt;
                    document.getElementById('es').style.color = 'red';
                }
            }
        };
    req.open("POST", es_script, true);
    req.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
    req.send(data);
    }

// ------ CHART FUNCTIONS -------
function update_chart(results){
    labels = [
            'PatmosO2',
            'PIO2',
            'PAO2',
            'PcO2',
            'PaO2',
            'PtO2',
            'PvO2',
        ]
    values = []
    for (index in labels){
        values[index] = results[labels[index]]
    }
    labels[0]='Atmos';//this looks neater
    cascadechart.data.labels = labels;
    cascadechart.data.datasets[0].data = values;
    cascadechart.update();

    var shift = 3.46877001647/results['P50_a']
    satschart.data.datasets[0].data = satline(shift, gasunits);
    satschart.update();

}

var ctx = document.getElementById('oxygencascadechart').getContext('2d');
var cascadechart = new Chart(ctx, { // https://www.chartjs.org/samples/latest/charts/line/basic.html
    type: 'bar',
    data: {
      labels: [
                "Atmos",
                "PIO2",
                "PAO2",
                "PcO2",
                "PaO2",
                "PtO2",
                "PvO2",
        ],
      datasets: [
        {
          label: "",
          backgroundColor: ["#3e95cd", "#3e95cd","#3e95cd","#3e95cd","#aa1212","#3e95cd","#3e95cd"],
          data: []
        }
      ]
    },
    options: {
        responsive: true,
        legend: { display: false },
        title: {
            display: false,
            text: 'O2cascade'
          },
        scales: {
            xAxes: [{
                barPercentage: 0.5,
                barThickness: 40,
                maxBarThickness: 60,
                minBarLength: 2,
                gridLines: {
                    offsetGridLines: true
                }
            }],
            yAxes: [{
                ticks: {
                    suggestedMin: 0,
                    suggestedMax: 30
                },
                scaleLabel: {
                    display: true,
                    labelString: 'Oxygen partial pressure (kPa)',
                    fontColor: 'gray',
                    padding: 0
                }
            }]
        }

    }
});