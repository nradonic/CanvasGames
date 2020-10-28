// CanvasGames - Annealing of colors - swapping two cell values if it is a better color match

var screenDraw = 0;
var screenDelay = 1;
var gridSize = 4;
var gridSize2 = gridSize * gridSize;
var serviceFlag = true; // set when operating parameters change
var Pause = true; // Pause Play button state
var Step = false; // runs through 1 cycle of smooth + FFT when true

var ColorSpace = 3;
var ColorScale = 255;
var maxTest = 3 * 255 * 255;

var avgPosn = {x:0, y:0};

//unique color area objects [{r,g,b}]
function ColorZone(red, green, blue) {
    this.red = red;
    this.green = green;
    this.blue = blue;
}

/**
 * Generate scaled color numbers - uses external color
 * @param colorNumber  maximum 1 less than colorSpace ie 0..1
 * @param colorSpace  number of points in the color palette   0...(n-1)   i.e. n=2
 * @returns {number}  final scaled number  ie 0...255 in 2 steps.... 0 and 255
 */
function scaleColor(colorNumber, colorSpace) {
    var tempColor = Math.min(colorNumber, colorSpace - 1);
    return Math.floor(tempColor * ColorScale / (colorSpace - 1));
}

// create array of color standards
function fillColorZones(colorSpace) {
    var colorZones = [];
    for (var r = 0; r < colorSpace; r++) {
        for (var g = 0; g < colorSpace; g++) {
            for (var b = 0; b < colorSpace; b++) {
                colorZones.push(new ColorZone(scaleColor(r, colorSpace), scaleColor(g, colorSpace), scaleColor(b, colorSpace)));
            }
        }
    }
    return colorZones;
}

// array of predetermined color objects
//unique color area objects [{r,g,b}]
var colorZones = fillColorZones(ColorSpace);

var changes = 0;
var totalChanges = 0;
var noiseLevel = 0;

var forceRange = 10;
var savedData = "";

var fftLayer = "All";

// var fixedSpots = [];
//
// function setFixedValues() {
//     // 4 x 3 array (row, column, linear_cell)
//     fixedSpots = [[Math.floor(gridSize / 4), Math.floor(gridSize / 4), Math.floor(gridSize / 4) * gridSize + Math.floor(gridSize / 4)],
//         [Math.floor(gridSize / 4), Math.floor(gridSize * 0.75), Math.floor(gridSize * 0.25) * gridSize + Math.floor(gridSize * 0.75)],
//         [Math.floor(gridSize * 0.75), Math.floor(gridSize / 4), Math.floor(gridSize * 0.75) * gridSize + Math.floor(gridSize * 0.25)],
//         [Math.floor(gridSize * 0.75), Math.floor(gridSize * 0.75), Math.floor(gridSize * 0.75) * gridSize + Math.floor(gridSize * 0.75)]];
// }
//
// setFixedValues();

// extension to Array type for 2D, with initialization - from Douglas Crockford
Array.matrix = function (numrows, numcols, initial) {
    var arr = [];
    for (var i = 0; i < numrows; ++i) {
        var columns = [];
        for (var j = 0; j < numcols; ++j) {
            columns[j] = {dist: 0, dist2: 0};
        }
        arr[i] = columns;
    }
    return arr;
};

function rand1() {
    return Math.floor(Math.random() * ColorSpace) * ColorScale;
}

function randGS(gridSizeX) {
    var a = Math.floor(Math.random() * gridSizeX);
    return a;
}

// calculate a +/- range variable
function randPMFr() {
    var a = Math.floor(Math.random() * (2 * forceRange + 1) - forceRange);
    return a;
}

function randPZM1() {   // returns -1, 0, +1
    var a = Math.floor(Math.random() * 3 - 1);
    return a;
}

function randX(x) {
    return math.floor(Math.random() * x);
}

function gridCell(rc, gc, bc) {
    this.r = rc;
    this.g = gc;
    this.b = bc;
    // keep indexes of colorZones....
    this.colorZoneIndex = randX(colorZones.size);
}

// dataGrid stores graphic data
var dataGrid = new Array(gridSize2);
var fftGrid = new Array(gridSize2);

// distanceGrid stores distance coefficients
var distanceGrid = Array.matrix(gridSize, gridSize, 0);

function fillDistanceGrid() {
    for (var i = 0; i < gridSize; i++) {
        for (var j = 0; j < gridSize; j++) {
            distanceGrid[i][j] = {dist: Math.sqrt(i * i + j * j), dist2: (i * i + j * j)};
        }
    }
    distanceGrid[0][0].dist = 1000;
    distanceGrid[0][0].dist2 = 1000000;
}

// fill dataGrid with color objects....
function fillDataGrid(dG, gs2) {
    for (var i = 0; i < gs2; i++) {
        dG[i] = new gridCell(rand1(), rand1(), rand1());
    }
}

// fill FFT grid with zeros in RGB elements
function zeroFFTGrid(fG, gS2) {
    for (var i = 0; i < gS2; i++) {
        fG[i] = new gridCell(0, 0, 0);
    }
}

// draw data
function drawData() {
    document.getElementById("screenDraw").innerHTML = screenDraw.toFixed(0);
    document.getElementById("changes").innerHTML = changes.toFixed(0);
    document.getElementById("totalChanges").innerHTML = totalChanges.toFixed(0);
    document.getElementById("noiseLevel").innerHTML = noiseLevel.toFixed(0);
}


// draw raw graphics pattern
function drawCanvas1() {
    var c = document.getElementById("drawHere");
    var ctx = c.getContext("2d");
    var myScreen = 800;
    ctx.beginPath();
    var squareSide = myScreen / gridSize;
    for (var i = 0; i < gridSize2; i++) {
        var squareRow = Math.floor(i / gridSize);
        var squareCol = Math.floor(i % gridSize);

        var y = squareRow * squareSide;
        var x = squareCol * squareSide;

        var t = 'rgba(' + dataGrid[i].colorZone.red + ',' + dataGrid[i].colorZone.green + ',' + dataGrid[i].colorZone.blue + ',255)';
        ctx.fillStyle = t;
        ctx.fillRect(x, y, squareSide, squareSide);
    }
    ctx.stroke();
}

// draw raw graphics pattern
function drawFFTCanvas() {
    var c = document.getElementById("drawFFT2D");
    var ctx = c.getContext("2d");
    var myScreen = 800;
    ctx.beginPath();
    var squareSide = myScreen / gridSize;
    for (var i = 0; i < gridSize2; i++) {
        var squareRow = Math.floor(i / gridSize);
        var squareCol = Math.floor(i % gridSize);

        var y = squareRow * squareSide;
        var x = squareCol * squareSide;
        var t = '';
        switch (fftLayer) {
            case "All":
                t = 'rgba(' + fftGrid[i].r + ',' + fftGrid[i].g + ',' + fftGrid[i].b + ',255)';
                break;
            case "Red":
                t = 'rgba(' + fftGrid[i].r + ',' + 0 + ',' + 0 + ',255)';
                break;
            case "Green":
                t = 'rgba(' + 0 + ',' + fftGrid[i].g + ',' + 0 + ',255)';
                break;
            case "Blue":
                t = 'rgba(' + 0 + ',' + 0 + ',' + fftGrid[i].b + ',255)';
                break;

            default:
                t = 'rgba(' + fftGrid[i].r + ',' + fftGrid[i].g + ',' + fftGrid[i].b + ',255)';
                break;
        }
        ctx.fillStyle = t;
        ctx.fillRect(x, y, squareSide, squareSide);
    }
    ctx.stroke();
}

function drawCanvas2(iIndex1, jIndex1, iIndex2, jIndex2) {
    var c = document.getElementById("drawHere");
    var ctx = c.getContext("2d");
    var myScreen = 800;
    ctx.beginPath();
    var squareSide = myScreen / gridSize;

    var i = iIndex1 * gridSize + jIndex1;
    var y = iIndex1 * squareSide;
    var x = jIndex1 * squareSide;


    var t = 'rgba(' + dataGrid[i].r + ',' + dataGrid[i].g + ',' + dataGrid[i].b + ',255)';
    ctx.fillStyle = t;
    ctx.fillRect(x, y, squareSide, squareSide);

    var j = iIndex2 * gridSize + jIndex2;
    var y2 = iIndex2 * squareSide;
    var x2 = jIndex2 * squareSide;
    t = 'rgba(' + dataGrid[j].r + ',' + dataGrid[j].g + ',' + dataGrid[j].b + ',255)';
    ctx.fillStyle = t;
    ctx.fillRect(x, y, squareSide, squareSide);

    ctx.lineWidth = 2;
    ctx.strokeStyle = '#FFFFFF';
    ctx.rect(x, y, squareSide, squareSide);
    ctx.rect(x2, y2, squareSide, squareSide);

    ctx.stroke();
}


// reset grid related values
function reconfigureGridRelatedStructures(dropdownSize, dropdownColor, dropdownRange) {
    gridSize = parseInt(dropdownSize.options[dropdownSize.selectedIndex].value);
    gridSize2 = gridSize * gridSize;
    dataGrid = new Array(gridSize2);
    fftGrid = new Array(gridSize2);

    // create new distance grid
    distanceGrid = Array.matrix(gridSize, gridSize, 0);
    // initialize distance grid
    fillDistanceGrid();

    // set the fixed color spot locations
    //setFixedValues();

    ColorSpace = parseInt(dropdownColor.options[dropdownColor.selectedIndex].value);
    ColorScale = Math.floor(255 / (ColorSpace - 1));

    cycleDelay();

    // generate new grid....
    fillDataGrid(dataGrid, gridSize2);
    zeroFFTGrid(fftGrid, gridSize2);
    screenDraw = 0;
}

// handle html button presses
function OnChange(param) {
    var dropdownSize = document.getElementById("select1");
    var dropdownColor = document.getElementById("select2");
    var dropdownRange = document.getElementById("select4");
    fftLayer = document.getElementById("selectFFTL").value;

    // adjust smoothing range
    if (param == 4) {
        forceRange = parseInt(dropdownRange.options[dropdownRange.selectedIndex].value);

    } else if (param == 5) {
        // already read in the new parameter fftLayer
        var a = 0; // for break point...
    } else {
        reconfigureGridRelatedStructures(dropdownSize, dropdownColor, dropdownRange);
    }
    drawData();
    drawCanvas1();
    fftGrid = fft2d(dataGrid, gridSize);
    drawFFTCanvas();
    serviceFlag = true;
}

function diffSQ0(ST1, ST2) {
    if (dataGrid[ST1].r === dataGrid[ST2].r &&
        dataGrid[ST1].g === dataGrid[ST2].g &&
        dataGrid[ST1].b === dataGrid[ST2].b) {
        return 0;
    } else {
        return 10000;
    }
}

// color space maximum less color differences squared
function diffSQ(ST1, ST2) {
    return maxTest - diffSQ0(ST1, ST2);
}

function fitColor(i1, j1, i2, j2) {
    var ST1 = i1 * gridSize + j1;
    var ST2 = i2 * gridSize + j2;
    var edge = gridSize - 1;

    var total = 0;
    var cnt = 0;
    if (j2 !== 0) {
        total = diffSQ(ST1, ST2 - 1);
        cnt++;
    }
    if (i2 !== 0 && j2 !== 0) {
        total += diffSQ(ST1, ST2 - gridSize - 1) / 2;
        cnt++;
    }//2;
    if (i2 !== 0) {
        total += diffSQ(ST1, ST2 - gridSize);
        cnt++;
    }
    if (i2 !== 0 && j2 !== edge) {
        total += diffSQ(ST1, ST2 - gridSize + 1) / 2;
        cnt++;
    }///2;
    if (j2 !== edge) {
        total += diffSQ(ST1, ST2 + 1);
        cnt++;
    }
    if (i2 !== edge && j2 !== 0) {
        total += diffSQ(ST1, ST2 + gridSize - 1) / 2;
        cnt++;
    }//2;
    if (i2 !== edge) {
        total += diffSQ(ST1, ST2 + gridSize);
        cnt++;
    }
    if (i2 !== edge && j2 !== edge) {
        total += diffSQ(ST1, ST2 + gridSize + 1) / 2;
        cnt++;
    }//2;
    return total / cnt;
}

//convert i position in array to column value x
function iToX(i, gS) {
    return Math.floor(i % gS);
}

//convert i position in array to row value y
function iToY(i, gS) {
    return Math.floor(i / gS);
}

// convert x,y to i in square grid
function xYToI(x, y, gS) {
    return x * gS + y;
}

// calculate weighted color-diff^2/distance^2 over entire array for a specific cell
function forceAtGridPoint(i1, j1) {
    var posn1 = xYToI(i1, j1, gridSize);
    // sum of forces variable
    var force = 0;

    var i1Min = Math.max(i1 - forceRange, 0);
    var i1Max = Math.min(i1 + forceRange, gridSize);

    var j1Min = Math.max(j1 - forceRange, 0);
    var j1Max = Math.min(j1 + forceRange, gridSize);

    for (var i = i1Min; i < i1Max; i++) {
        for (var j = j1Min; j < j1Max; j++) {
            // get grid locations for 'i' element   2D->1D
            // calculate position differences
            var dx = j1 - j;
            var dy = i1 - i;
            // lookup location difference squared vector
            var dist2 = distanceGrid[Math.abs(dy)][Math.abs(dx)].dist2;
            var posn0 = xYToI(i, j, gridSize);
            // calculate color difference squared
            var dc = diffSQ0(posn0, posn1);
            // add weighted force to sum of forces
            var forceScalar = dc * dist2;
            force += forceScalar;
        }
    }
    return force;
}

// trade colors at two locations
function trade(i1, j1, i2, j2) {
    var one = i1 * gridSize + j1;
    var two = i2 * gridSize + j2;

    var r = dataGrid[one].r;
    var g = dataGrid[one].g;
    var b = dataGrid[one].b;

    dataGrid[one].r = dataGrid[two].r;
    dataGrid[one].g = dataGrid[two].g;
    dataGrid[one].b = dataGrid[two].b;

    dataGrid[two].r = r;
    dataGrid[two].g = g;
    dataGrid[two].b = b;
}

// returns force reduction for swap : positive is an improvement, negative is a worsening of differences
function fitColorAll(i1, j1, i2, j2, force1) {
    var posn1 = xYToI(i1, j1, gridSize);
    var posn2 = xYToI(i2, j2, gridSize);
    //var force1 = forceAtGridPoint(i1,j1);
    var force2 = forceAtGridPoint(i2, j2);
    trade(i1, j1, i2, j2);
    var force1At2 = forceAtGridPoint(i2, j2);
    var force2At1 = forceAtGridPoint(i1, j1);
    trade(i1, j1, i2, j2);
    return force1 + force2 - force1At2 - force2At1;
}

// Test swaps with all nearest neighbors
function testCellNeighborsOverGrid(i1, j1) {
    var eta = 1;

    // force plus 1 row plus one column
    function fPP(i1, j1, force1At1) {
        if (i1 + 1 < gridSize && j1 + 1 < gridSize) {
            return fitColorAll(i1, j1, i1 + 1, j1 + 1, force1At1);
        } else {
            return 0;
        }
    }

    // force plus 1 row minus one column
    function fPN(i1, j1, force1At1) {
        if (i1 + 1 < gridSize && j1 - 1 >= 0) {
            return fitColorAll(i1, j1, i1 + 1, j1 - 1, force1At1);
        } else {
            return 0;
        }
    }

    // force plus 1 row same column
    function fPZ(i1, j1, force1At1) {
        if (i1 + 1 < gridSize) {
            return fitColorAll(i1, j1, i1 + 1, j1, force1At1);
        } else {
            return 0;
        }
    }

    // force minus 1 row plus one column
    function fNP(i1, j1, force1At1) {
        if (i1 - 1 >= 0 && j1 + 1 < gridSize) {
            return fitColorAll(i1, j1, i1 - 1, j1 + 1, force1At1);
        } else {
            return 0;
        }
    }

    // force minus 1 row minus one column
    function fNN(i1, j1, force1At1) {
        if (i1 - 1 >= 0 && j1 - 1 >= 0) {
            return fitColorAll(i1, j1, i1 - 1, j1 - 1, force1At1);
        } else {
            return 0;
        }
    }

    // force minus 1 row same column
    function fNZ(i1, j1, force1At1) {
        if (i1 - 1 >= 0) {
            return fitColorAll(i1, j1, i1 - 1, j1, force1At1);
        } else {
            return 0;
        }
    }

    // force same 1 row minus one column
    function fZN(i1, j1, force1At1) {
        if (j1 - 1 >= 0) {
            return fitColorAll(i1, j1, i1, j1 - 1, force1At1);
        } else {
            return 0;
        }
    }

    // force same row plus one column
    function fZP(i1, j1, force1At1) {
        if (j1 + 1 < gridSize) {
            return fitColorAll(i1, j1, i1, j1 + 1, force1At1);
        } else {
            return 0;
        }
    }

    var force1At1 = forceAtGridPoint(i1, j1);
    // as the 'force' increases there is more advantage to swapping...
    var maxForce = 1;  // minimum noise margin
    var minVect = "ZZ";
    var ffPP = fPP(i1, j1, force1At1);
    if (ffPP > maxForce) {
        maxForce = ffPP;
        minVect = "PP";
    }
    var ffPN = fPN(i1, j1, force1At1);
    if (ffPN > maxForce) {
        maxForce = ffPN;
        minVect = "PN";
    }
    var ffPZ = fPZ(i1, j1, force1At1);
    if (ffPZ > maxForce) {
        maxForce = ffPZ;
        minVect = "PZ";
    }
    var ffNN = fNN(i1, j1, force1At1);
    if (ffNN > maxForce) {
        maxForce = ffNN;
        minVect = "NN";
    }
    var ffNP = fNP(i1, j1, force1At1);
    if (ffNP > maxForce) {
        maxForce = ffNP;
        minVect = "NP";
    }
    var ffNZ = fNZ(i1, j1, force1At1);
    if (ffNZ > maxForce) {
        maxForce = ffNZ;
        minVect = "NZ";
    }
    var ffZP = fZP(i1, j1, force1At1);
    if (ffZP > maxForce) {
        maxForce = ffZP;
        minVect = "ZP";
    }
    var ffZN = fZN(i1, j1, force1At1);
    if (ffZN > maxForce) {
        maxForce = ffZN;
        minVect = "ZN";
    }

    if (minVect !== "ZZ") {
        var lchanges = changes;
        changes++;
        if (maxForce < 0) {
            var m = maxForce;
        }
        noiseLevel += maxForce;
    }

    switch (minVect) {
        case "PP":
            trade(i1, j1, i1 + 1, j1 + 1);
            break;
        case "PN":
            trade(i1, j1, i1 + 1, j1 - 1);
            break;
        case "PZ":
            trade(i1, j1, i1 + 1, j1);
            break;
        case "NN":
            trade(i1, j1, i1 - 1, j1 - 1);
            break;
        case "NP":
            trade(i1, j1, i1 - 1, j1 + 1);
            break;
        case "NZ":
            trade(i1, j1, i1 - 1, j1);
            break;
        case "ZP":
            trade(i1, j1, i1, j1 + 1);
            break;
        case "ZN":
            trade(i1, j1, i1, j1 - 1);
            break;
        default: {
            break;
        }
    }
}

// scan for swap over entire grid space with vector force measurements
function vectorSwap() {
    changes = 0;
    noiseLevel = 0;
    var eta = 1;

    for (var i = 0; i < gridSize2; i++) {
        var i1 = iToY(i, gridSize);
        var j1 = iToX(i, gridSize);
        lchanges = changes;
        testCellNeighborsOverGrid(i1, j1);
    }
    // Annealing - random cell swaps....
    var gSM1 = gridSize - 1;
    for (var i = 0; i < gridSize; i++) {
        var i1 = randGS(gridSize);   // i,j +/- foreceRange....
        var j1 = randGS(gridSize);
        var force1At1 = forceAtGridPoint(i1, j1);
        var i2 = Math.min(gSM1, Math.max(0, i1 + randPMFr()));
        var j2 = Math.min(gSM1, Math.max(0, j1 + randPMFr()));
        var force2At2 = forceAtGridPoint(i2, j2);
        trade(i1, j1, i2, j2);
        var force1At2 = forceAtGridPoint(i2, j2);
        var force2At1 = forceAtGridPoint(i1, j1);
        var resetTrade = force1At1 + force2At2 + eta - force1At2 - force2At1;
        if (resetTrade < 0) {
            trade(i1, j1, i2, j2);  // reset if no gain in trades...
        } else {
            lchanges = changes;
            changes++;
            noiseLevel += resetTrade;
        }
    }
    lchanges = changes;
    totalChanges += changes;
}

var cellColors = {
    colors: [],
    findOrAddColor: function (cellColor, cell, i, j) {
        var found = false;
        this.colors.forEach(function (element) {
            if (element.r === cellColor.r && element.g === cellColor.g && element.b === cellColor.b) {
                found = true;
            }
        })
        if (!found) {
            this.colors.push(cellColor);
        }
    }
};

var steps = [1, -1, gridSize, -gridSize];

var colorCenters = [];

/**
 * zeroColorCenters - zeros out the x,y location parameters for all color zones
 * @param colorZones
 */
function zeroColorCenters(colorZones){
    colorCenters = [];
    for (var i=0; i<colorZones.size; i++){
        colorCenters.push({x:0, y:0});
    }
}

/**
 * calculates the average x and y centroid for each color pattern
 *
 */
function calculateColorCenters(){
    // fill the colors array - one object per colorZone pattern
    colorCenters = zeroColorCenters(colorZones);
    //retrieve color index per cell, increment those position sums
    for (var x=0; x< gridSize; x++){
        for ( var y=0; y<gridSize; y++){
            var pos = x * gridSize + y;
            // get color index at each grid point
            var colorGridPt = dataGrid[pos].colorZone;

            var z = colorCenters[colorGridPt];
            z.x += x;
            z.y += y;
            colorCenters[colorGridPt] = z;
        }
    }
}

function vectorSwap2() {
    changes = 0;
    noiseLevel = 0;
    var eta = 1;

    calculateColorCenters();

    var lchanges = changes;
    for (var i = 0; i < gridSize2; i++) {
        var i1 = iToY(i, gridSize);
        var j1 = iToX(i, gridSize);
        lchanges = changes;
        testCellNeighborsOverGrid(i1, j1);
    }
    // Annealing - random cell swaps....
    var gSM1 = gridSize - 1;
    for (var i = 0; i < gridSize; i++) {
        var i1 = randGS(gridSize);
        var j1 = randGS(gridSize);
        var force1At1 = forceAtGridPoint(i1, j1);

        var i2 = Math.max(Math.min(i1 + step[randX(2)], gridSize), 0);
        var j2 = Math.max(0, Math.min(j1 + steps[randX(2)], gridSize));
        var force2At2 = forceAtGridPoint(i2, j2);
        trade(i1, j1, i2, j2);
        var force1At2 = forceAtGridPoint(i2, j2);
        var force2At1 = forceAtGridPoint(i1, j1);
        var resetTrade = force1At1 + force2At2 + eta - force1At2 - force2At1;
        if (resetTrade < 0) {
            trade(i1, j1, i2, j2);  // reset if no gain in trades...
        } else {
            lchanges = changes;
            changes++;
            noiseLevel += resetTrade;
        }
    }
    totalChanges += changes;
}

function smooth() {
    screenDraw++;
    drawData();

    vectorSwap2();
    //vectorSwap2();
    fftGrid = fft2d(dataGrid, gridSize);
}

function cycle(myVarr) {
    smooth();
    drawCanvas1();
    drawFFTCanvas();
    if (serviceFlag) {
        clearInterval(myVarr);
        if (Step === false) {
            startLooping();
        }
    }
}

var myVar;

// start calculating updates
function startLooping() {
    myVar = setInterval(function () {
        cycle(myVar)
    }, screenDelay);
    Pause = false;
    document.getElementById("PausePlay").innerHTML = "Pause";
}

// start repeating operations
function start() {
    Step = false;
    startLooping();
    //document.getElementById("PausePlay").innerHTML="Pause";
}

function stop() {
    clearInterval(myVar);
    Pause = true;
    document.getElementById("PausePlay").innerHTML = "Play";
}

function step() {
    stop();
    Step = true;
    cycle(myVar);
}

function cycleDelay() {
    var dropdownDelay = document.getElementById("select3");
    screenDelay = parseInt(dropdownDelay.options[dropdownDelay.selectedIndex].value);
    serviceFlag = true;
}

function ButtonLabelToPlay() {
    document.getElementById("PausePlay").innerHTML = "Play";
}

// toggle between running and paused if current state is not paused
function PausePlay() {
    if (Pause) {
        start();
    }
    else {
        stop();
    }
}

//fillDistanceGrid();
OnChange(0);
//drawCanvas1();
//start();

function SaveData() {
    stop();
    savedData = JSON.stringify({
        gridSize: gridSize,
        gridSize2: gridSize2,
        screenDraw: screenDraw,
        ColorSpace: ColorSpace,
        screenDelay: screenDelay,
        dataGrid: dataGrid
    });
    ButtonLabelToPlay();

    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/Controller/Action");
    xhr.setRequestHeader('Content-Type', 'application/json');
    xhr.onreadystatechange = function () {
        if (xhr.readyState == 4 && xhr.status == 200) {
            alert(xhr.responseText);
        }
    }
    xhr.send(savedData);
}

function LoadData() {
    stop();
    var b = JSON.parse(savedData);
    gridSize = b.gridSize;
    gridSize2 = b.gridSize2;
    screenDraw = b.screenDraw;
    screenDelay = b.screenDelay;
    ColorSpace = b.ColorSpace;
    dataGrid = b.dataGrid;
    drawCanvas1();
    ButtonLabelToPlay();
}