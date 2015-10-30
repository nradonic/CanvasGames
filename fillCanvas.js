// CanvasGames - Annealing of colors - swapping two cell values if it is a better color match

var screenDraw = 0;
var screenDelay=1;
var gridSize = 4;
var gridSize2 = gridSize*gridSize;
var serviceFlag = true; // set when operating parameters change
var Pause = true; // Pause Play button state
var Step = false; // runs through 1 cycle of smooth + FFT when true

var ColorSpace = 2;
var ColorScale = 255;
var maxTest = 3*255*255;
var screenDrawI = -1;
var screenDrawJ = -1;

var changes = 0;
var totalChanges = 0;
var noiseLevel = 0;

var forceRange = 10;
var savedData="";

var dbgMinST1 = 10000000;
var dbgMaxST1 = -10000000;
var dbgMinST2 = 10000000;
var dbgMaxST2 = -10000000;

var fftLayer = "All";

var fixedSpots = [];

function setFixedValues(){
	// 4 x 3 array (row, column, linear_cell)
	fixedSpots = [[Math.floor(gridSize/4), Math.floor(gridSize/4), Math.floor(gridSize/4)*gridSize+Math.floor(gridSize/4)],
					[Math.floor(gridSize/4), Math.floor(gridSize*0.75), Math.floor(gridSize*0.25)*gridSize+Math.floor(gridSize*0.75)],
					[Math.floor(gridSize*0.75), Math.floor(gridSize/4), Math.floor(gridSize*0.75)*gridSize+Math.floor(gridSize*0.25)],
					[Math.floor(gridSize*0.75), Math.floor(gridSize*0.75), Math.floor(gridSize*0.75)*gridSize+Math.floor(gridSize*0.75)]];	
}

setFixedValues();

// extension to Array type for 2D, with initialization - from Douglas Crockford
Array.matrix = function(numrows, numcols, initial){
   var arr = [];
   for (var i = 0; i < numrows; ++i){
      var columns = [];
      for (var j = 0; j < numcols; ++j){
         columns[j] = {dist:0, dist2:0};
      }
      arr[i] = columns;
    }
    return arr;
};

function rand1(){
	return Math.floor(Math.random()*ColorSpace)*ColorScale;
}

function randGS(gridSizeX){
	var a = Math.floor(Math.random()*gridSizeX);
	return a;
}

// calculate a +/- range variable
function randPMFr(){ 
	var a = Math.floor(Math.random()*(2*forceRange+1)-forceRange);
	return a;
}

function randPZM1(){   // returns -1, 0, +1
	var a = Math.floor(Math.random()*3 -1);
	return a;
}

function gridCell(rc,gc,bc){this.r=rc;this.g=gc;this.b=bc;}

// dataGrid stores graphic data
var dataGrid = new Array(gridSize2);
var fftGrid = new Array(gridSize2);

// distanceGrid stores distance coefficients
var distanceGrid = Array.matrix(gridSize,gridSize,0);

function fillDistanceGrid(){
	for (var i=0; i<gridSize; i++){
		for (var j=0; j<gridSize; j++){
			distanceGrid[i][j] = {dist:Math.sqrt(i*i + j*j), dist2:(i*i + j*j)};
		}
	}
	distanceGrid[0][0].dist=1000;
	distanceGrid[0][0].dist2=1000000;
}

// fill dataGrid with color objects....
function fillDataGrid(dG, gs2){
	for(var i=0;i<gs2;i++){
//  		dG[i]=new gridCell(0,0,0);
	//if(i%3===0){
	//	dG[i]=new gridCell(0,0,0);
	//} else {
		dG[i]=new gridCell(rand1(),rand1(),rand1());
	//	dG[i]=new gridCell(255,0,0);
	//} 

	}
//	dG[fixedSpots[0][2]] = new gridCell(255,0,0); // red
//	dG[fixedSpots[1][2]] = new gridCell(0,255,0); // green
//	dG[fixedSpots[2][2]] = new gridCell(0,0,255); // blue
//	dG[fixedSpots[3][2]] = new gridCell(0,0,0);   // black
}

// fill FFT grid with zeros in RGB elements
function zeroFFTGrid(fG, gS2){
	for ( var i = 0; i < gS2; i++){
		fG[i] = new gridCell(0,0,0);
	}
}


// draw raw graphics pattern
function drawCanvas1(){
	var c=document.getElementById("drawHere");
	var ctx=c.getContext("2d");
	var myScreen = 800;
	var myScreen2 = myScreen*myScreen*4;
	var myScreen4 = myScreen*4;
	var scaler =myScreen/gridSize; 
	ctx.beginPath();
	var squareSide = myScreen/gridSize;
	for(var i = 0; i< gridSize2;i++){
		var squareRow = Math.floor(i/gridSize);
		var squareCol = Math.floor(i%gridSize);

		var y = squareRow*squareSide;
		var x = squareCol*squareSide;
		
		var t='rgba('+dataGrid[i].r+','+ dataGrid[i].g+','+dataGrid[i].b+',255)';
		ctx.fillStyle = t;
		ctx.fillRect(x,y,squareSide,squareSide);
    }
   	ctx.stroke();
}


// draw raw graphics pattern
function drawFFTCanvas(){
	var c=document.getElementById("drawFFT2D");
	var ctx=c.getContext("2d");
	var myScreen = 800;
	var myScreen2 = myScreen*myScreen*4;
	var myScreen4 = myScreen*4;
	var scaler =myScreen/gridSize; 
	ctx.beginPath();
	var squareSide = myScreen/gridSize;
	for(var i = 0; i< gridSize2;i++){
		var squareRow = Math.floor(i/gridSize);
		var squareCol = Math.floor(i%gridSize);

		var y = squareRow*squareSide;
		var x = squareCol*squareSide;
		var t= '';
		switch(fftLayer){
			case "All":
				t='rgba('+fftGrid[i].r+','+ fftGrid[i].g+','+fftGrid[i].b+',255)';
				break;		
			case "Red":
				t='rgba('+fftGrid[i].r+','+ 0+','+0+',255)';
				break;		
			case "Green":
				t='rgba('+0+','+ fftGrid[i].g+','+0+',255)';
				break;		
			case "Blue":
				t='rgba('+0+','+ 0+','+fftGrid[i].b+',255)';
				break;		
			 
		    default:
  				t='rgba('+fftGrid[i].r+','+ fftGrid[i].g+','+fftGrid[i].b+',255)';
		    	break;
		}
		ctx.fillStyle = t;
		ctx.fillRect(x,y,squareSide,squareSide);
    }
   	ctx.stroke();
}

function drawCanvas2(iIndex1, jIndex1, iIndex2, jIndex2){
	var c=document.getElementById("drawHere");
	var ctx=c.getContext("2d");
	var myScreen = 800;
	var myScreen2 = myScreen*myScreen*4;
	var myScreen4 = myScreen*4;
	var scaler =myScreen/gridSize; 
	ctx.beginPath();
	var squareSide = myScreen/gridSize;
	
	var i = iIndex1*gridSize + jIndex1;
	var y = iIndex1*squareSide;
	var x = jIndex1*squareSide;
		
		
	var t='rgba('+dataGrid[i].r+','+ dataGrid[i].g+','+dataGrid[i].b+',255)';
	ctx.fillStyle = t;
	ctx.fillRect(x,y,squareSide,squareSide);
	
	var j = iIndex2*gridSize + jIndex2;
	var y2 = iIndex2*squareSide;
	var x2 = jIndex2*squareSide;
	t='rgba('+dataGrid[j].r+','+ dataGrid[j].g+','+dataGrid[j].b+',255)';
	ctx.fillStyle = t;
	ctx.fillRect(x,y,squareSide,squareSide);

    ctx.lineWidth=2;
    ctx.strokeStyle= '#FFFFFF';
    ctx.rect(x,y,squareSide,squareSide);
    ctx.rect(x2,y2,squareSide,squareSide);
    
	ctx.stroke();
}


// reset grid related values
function reconfigureGridRelatedStructures(dropdownSize, dropdownColor, dropdownRange){
	    gridSize = parseInt(dropdownSize.options[dropdownSize.selectedIndex].value);
	    gridSize2 = gridSize*gridSize;   
	    dataGrid = new Array(gridSize2);
	    fftGrid = new Array(gridSize2);
	    
	    // create new distance grid
	    distanceGrid = Array.matrix(gridSize,gridSize,0);
	    // initialize distance grid 
	    fillDistanceGrid();

	    // set the fixed color spot locations
	    setFixedValues();

	    ColorSpace = parseInt(dropdownColor.options[dropdownColor.selectedIndex].value);
	    ColorScale = Math.floor(255/(ColorSpace-1));
	    
	    cycleDelay();

        // generate new grid....
	    fillDataGrid(dataGrid, gridSize2);	    
	    zeroFFTGrid(fftGrid, gridSize2);
        screenDraw = 0;
}

// handle html button presses
function OnChange(param)
{
	var dropdownSize = document.getElementById("select1");
    var dropdownColor = document.getElementById("select2");
    var dropdownRange = document.getElementById("select4");
    fftLayer = document.getElementById("selectFFTL").value;

    // adjust smoothing range
    if(param==4){
    	forceRange = parseInt(dropdownRange.options[dropdownRange.selectedIndex].value);
    
    } else if(param==5){
    	// already read in the new parameter fftLayer
    	var a=0; // for break point...
    } else {
// 	    gridSize = parseInt(dropdownSize.options[dropdownSize.selectedIndex].value);
// 	    gridSize2 = gridSize*gridSize;   
// 	    dataGrid = new Array(gridSize2);
// 	    fftGrid = new Array(gridSize2);
// 	    
// 	    // create new distance grid
// 	    distanceGrid = Array.matrix(gridSize,gridSize,0);
// 	    // initialize distance grid 
// 	    fillDistanceGrid();
// 
// 	    // set the fixed color spot locations
// 	    setFixedValues();
// 
// 	    ColorSpace = parseInt(dropdownColor.options[dropdownColor.selectedIndex].value);
// 	    ColorScale = Math.floor(255/(ColorSpace-1));
// 	    
// 	    cycleDelay();
// 
//         // generate new grid....
// 	    fillDataGrid(dataGrid, gridSize2);	    
// 	    zeroFFTGrid(fftGrid, gridSize2);
//         screenDraw = 0;

		reconfigureGridRelatedStructures(dropdownSize, dropdownColor, dropdownRange);
    }
	drawCanvas1(); 
	drawFFTCanvas();
    serviceFlag = true;
    
//    return true;
}


// color differences squared only

function diffSQ0(ST1,ST2){
		dbgMinST1 = Math.min(dbgMinST1,ST1);
		dbgMinST2 = Math.min(dbgMinST2,ST2);
		dbgMaxST1 = Math.max(dbgMaxST1,ST1);
		dbgMaxST2 = Math.max(dbgMaxST2,ST2);
		var rd = dataGrid[ST1].r-dataGrid[ST2].r;
		var gd = dataGrid[ST1].g-dataGrid[ST2].g;
		var bd = dataGrid[ST1].b-dataGrid[ST2].b;
		return rd*rd+gd*gd+bd*bd;
}

// color space maximum less color differences squared
function diffSQ(ST1,ST2){
		return maxTest - diffSQ0(ST1,ST2);
}

function fitColor(i1,j1,i2,j2){
	var ST1 = i1*gridSize+j1;
	var ST2 = i2*gridSize+j2;
	var edge = gridSize-1;
	
	var total = 0;
	var cnt = 0;
	if(j2!== 0){total = diffSQ(ST1, ST2-1); cnt++;}
	if(i2!==0 && j2!==0){total += diffSQ(ST1,ST2-gridSize-1)/2; cnt++;}//2;
	if(i2!==0){total += diffSQ(ST1,ST2-gridSize); cnt++;}
	if(i2!==0 && j2!==edge){total += diffSQ(ST1,ST2-gridSize+1)/2; cnt++;}///2;
	if(j2!==edge){total += diffSQ(ST1,ST2+1); cnt++;}
	if(i2!==edge && j2!==0){total += diffSQ(ST1,ST2+gridSize-1)/2; cnt++;}//2;
	if(i2!==edge){total += diffSQ(ST1,ST2+gridSize); cnt++;}
	if(i2!==edge && j2!==edge){total += diffSQ(ST1,ST2+gridSize+1)/2; cnt++;}//2;
	return total/cnt;
}

//convert i position in array to column value x
function iToX(i,gS){
	return Math.floor(i%gS);
}

//convert i position in array to row value y
function iToY(i,gS){
	return Math.floor(i/gS);
}
// convert x,y to i in square grid
function xYToI(x,y,gS){
	return x*gS+y;	
}

// calculate weighted color-diff^2/distance^2 over entire array for a specific cell
function forceAtGridPoint(i1,j1){
	var posn1 = xYToI(i1,j1,gridSize);
		// sum of forces variable
	var force = 0;
	
	var i1Min= Math.max(i1-forceRange,0);
	var i1Max= Math.min(i1+forceRange,gridSize);
	
	var j1Min= Math.max(j1-forceRange,0);
	var j1Max= Math.min(j1+forceRange,gridSize);
	
	for (var i= i1Min; i<i1Max; i++){
		for (var j= j1Min; j<j1Max; j++){
			// get grid locations for 'i' element   2D->1D
		    // calculate position differences
			var dx = j1-j;
			var dy = i1-i;
			// lookup location difference squared vector
			var dist2 = distanceGrid[Math.abs(dy)][Math.abs(dx)].dist2;
			var posn0 = xYToI(i,j,gridSize);
			// calculate color difference squared
			var dc = diffSQ0(posn0, posn1);
			// add weighted force to sum of forces
			var forceScalar = dc / dist2;
			force += forceScalar;
		}
	}
	return force;
}

// trade colors at two locations
function trade(i1,j1,i2,j2){
	var one = i1*gridSize+j1;
	var two = i2*gridSize+j2;
	
	var r = dataGrid[one].r;
	var g = dataGrid[one].g;
	var b = dataGrid[one].b;
	
	dataGrid[one].r=dataGrid[two].r;
	dataGrid[one].g=dataGrid[two].g;
	dataGrid[one].b=dataGrid[two].b;
	
	dataGrid[two].r=r;
	dataGrid[two].g=g;
	dataGrid[two].b=b;
}	

// returns force reduction for swap : positive is an improvement, negative is a worsening of differences
function fitColorAll(i1,j1,i2,j2,force1){
	var posn1 = xYToI(i1,j1,gridSize);
	var posn2 = xYToI(i2,j2,gridSize);
	//var force1 = forceAtGridPoint(i1,j1);
	var force2 = forceAtGridPoint(i2,j2);
	trade(i1,j1,i2,j2);
	var force1At2 = forceAtGridPoint(i2,j2);
	var force2At1 = forceAtGridPoint(i1,j1);
	trade(i1,j1,i2,j2);
	return force1 + force2 - force1At2 - force2At1;
}

// Test swaps with all nearest neighbors
function testCellNeighborsOverGrid(i1,j1){
	var eta = 1;
	// force plus 1 row plus one column
	function fPP(i1,j1,force1At1){
		if(i1+1<gridSize && j1+1<gridSize){
			return fitColorAll(i1,j1,i1+1,j1+1,force1At1);
		} else {return 0;}
	}
	
	// force plus 1 row minus one column
	function fPN(i1,j1,force1At1){
		if(i1+1<gridSize && j1-1>=0){
			return fitColorAll(i1,j1,i1+1,j1-1,force1At1);
		} else {return 0;}
	}
	
	// force plus 1 row same column
	function fPZ(i1,j1,force1At1){
		if(i1+1<gridSize){
			return fitColorAll(i1,j1,i1+1,j1,force1At1);
		} else {return 0;}
	}
	
	// force minus 1 row plus one column
	function fNP(i1,j1,force1At1){
		if(i1-1>=0 && j1+1<gridSize){
			return fitColorAll(i1,j1,i1-1,j1+1,force1At1);
		} else {return 0;}
	}
	
	// force minus 1 row minus one column
	function fNN(i1,j1,force1At1){
		if(i1-1>=0 && j1-1>=0){
			return fitColorAll(i1,j1,i1-1,j1-1,force1At1);
		} else {return 0;}
	}
	
	// force minus 1 row same column
	function fNZ(i1,j1,force1At1){
		if(i1-1>=0){
			return fitColorAll(i1,j1,i1-1,j1,force1At1);
		} else {return 0;}
	}
	
	// force same 1 row minus one column
	function fZN(i1,j1,force1At1){
		if(j1-1>=0){
			return fitColorAll(i1,j1,i1,j1-1,force1At1);
		} else {return 0;}
	}
	
	// force same row plus one column
	function fZP(i1,j1,force1At1){
		if(j1+1<gridSize){
			return fitColorAll(i1,j1,i1,j1+1,force1At1);
		} else {return 0;}
	}

	var force1At1 = forceAtGridPoint(i1,j1);
	// as the 'force' increases there is more advantage to swapping...
	var maxForce = 1;  // minimum noise margin
	var minVect = "ZZ";
	var ffPP = fPP(i1,j1,force1At1);
	if(ffPP>maxForce){
		maxForce = ffPP;
		minVect = "PP";
	}
	var ffPN = fPN(i1,j1,force1At1);
	if(ffPN>maxForce){
		maxForce = ffPN;
		minVect = "PN";
	}
	var ffPZ = fPZ(i1,j1,force1At1);
	if(ffPZ>maxForce){
		maxForce = ffPZ;
		minVect = "PZ";
	}	
	var ffNN = fNN(i1,j1,force1At1);
	if(ffNN>maxForce){
		maxForce = ffNN;
		minVect = "NN";
	}
	var ffNP = fNP(i1,j1,force1At1);
	if(ffNP>maxForce){
		maxForce = ffNP;
		minVect = "NP";
	}
	var ffNZ = fNZ(i1,j1,force1At1);
	if(ffNZ>maxForce){
		maxForce = ffNZ;
		minVect = "NZ";
	}
	var ffZP = fZP(i1,j1,force1At1);
	if(ffZP>maxForce){
		maxForce = ffZP;
		minVect = "ZP";
	}
	var ffZN = fZN(i1,j1,force1At1);
	if(ffZN>maxForce){
		maxForce = ffZN;
		minVect = "ZN";
	}
	
	if(minVect!=="ZZ") { 
		var lchanges = changes;
		changes++;
		if (maxForce < 0){ 
		 var m = maxForce;
		}
		noiseLevel += maxForce;
	}
	
	switch(minVect){
		case "PP":
			trade(i1,j1,i1+1,j1+1);
			break;
		case "PN":
			trade(i1,j1,i1+1,j1-1);
			break;
		case "PZ":
			trade(i1,j1,i1+1,j1);
			break;
		case "NN":
			trade(i1,j1,i1-1,j1-1);
			break;
		case "NP":
			trade(i1,j1,i1-1,j1+1);
			break;
		case "NZ":
			trade(i1,j1,i1-1,j1);
			break;
		case "ZP":
			trade(i1,j1,i1,j1+1);
			break;
		case "ZN":
			trade(i1,j1,i1,j1-1);
			break;
	    default:
	    	{
	    		break;
	    	}
	}
}
	
// scan for swap over entire grid space with vector force measurements
function vectorSwap(){
	changes = 0;
	noiseLevel = 0;
	var eta = 1;
	
	var lchanges = changes;
	var a = fixedSpots[0][2];
	var b = fixedSpots[1][2];
	var c = fixedSpots[2][2];
	var d = fixedSpots[3][2];
	
	for (var i=0;i<gridSize2;i++){
		var i1 = iToY(i,gridSize);
		var j1 = iToX(i,gridSize);
		lchanges = changes;
			testCellNeighborsOverGrid(i1,j1);
	}
	// Annealing - random cell swaps....
	var gSM1 = gridSize-1;
	for (var i = 0; i<gridSize; i++){
		var i1 = randGS(gridSize);   // i,j +/- foreceRange....
		var j1 = randGS(gridSize);
		var force1At1 = forceAtGridPoint(i1,j1);
		var i2 = Math.min(gSM1,Math.max(0,i1 + randPMFr()));
		var j2 = Math.min(gSM1,Math.max(0,j1 + randPMFr()));
		var force2At2 = forceAtGridPoint(i2,j2);
		trade(i1,j1,i2,j2);
		var force1At2 = forceAtGridPoint(i2,j2);
		var force2At1 = forceAtGridPoint(i1,j1);
		var resetTrade = force1At1 + force2At2 + eta - force1At2 - force2At1;
		if(resetTrade < 0){
			trade(i1,j1,i2,j2);  // reset if no gain in trades...
		} else {
			lchanges = changes;
			changes++;
			noiseLevel += resetTrade;
		}
	}
	lchanges = changes;
	totalChanges += changes;
}	
	
function testSwapSpaces(i1,j1,i2,j2){
	var OneAtOne = fitColor(i1,j1,i1,j1);
	var TwoAtTwo = fitColor(i2,j2,i2,j2);
	var OneAtTwo = fitColor(i1,j1,i2,j2);
	var TwoAtOne = fitColor(i2,j2,i1,j1);
	if(Math.max(OneAtOne,TwoAtTwo)<Math.max(OneAtTwo,TwoAtOne)){trade(i1,j1,i2,j2);}
}	

function middleSwap(){
	for (var i=0;i<gridSize2;i++){
		setTimeout(function(){},screenDelay);
		var iIndex1 = Math.floor(i/gridSize); //randGS(gridSize);
		var jIndex1 = Math.floor(i%gridSize); //randGS(gridSize);
		var iIndex2 = Math.max(0,Math.min(gridSize-1,iIndex1 + randPZM1()));
		var jIndex2 = Math.max(0,Math.min(gridSize-1,jIndex1 + randPZM1()));
		if(iIndex1!==iIndex2 && jIndex1!==jIndex2){
			testSwapSpaces(iIndex1,jIndex1,iIndex2,jIndex2);
			drawCanvas2(iIndex1, jIndex1, iIndex2, jIndex2);
		}
	 }
}

function smooth(){
	screenDraw++;
	document.getElementById("screenDraw").innerHTML=screenDraw.toFixed(0);
	document.getElementById("changes").innerHTML=changes.toFixed(0);
	document.getElementById("totalChanges").innerHTML=totalChanges.toFixed(0);
	document.getElementById("noiseLevel").innerHTML=noiseLevel.toFixed(0);
	
	//middleSwap();
	vectorSwap();
	fftGrid = fft2d(dataGrid, gridSize);
}

function cycle(myVarr){
	drawCanvas1(); 
	drawFFTCanvas();
	smooth();
	if(serviceFlag){
		clearInterval(myVarr);
		if(Step===false){
				startLooping();
		}
	}
}

var myVar;


function startLooping(){
	myVar=setInterval(function(){cycle(myVar)},screenDelay);
	Pause = false;
}

// start repeating operations
function start(){
	Step = false;
	startLooping();
}

function stop(){
	clearInterval(myVar); 
	Pause = true;
}

function step(){
	Step = true;
	cycle(myVar);
}

function cycleDelay(){
	var dropdownDelay = document.getElementById("select3");
	screenDelay = parseInt(dropdownDelay.options[dropdownDelay.selectedIndex].value);
	serviceFlag = true;
}

function ButtonLabelToPlay(){
	document.getElementById("PausePlay").innerHTML="Play";
}

function PausePlay(){
	Pause = !Pause;
	if(Pause){stop(); document.getElementById("PausePlay").innerHTML="Play";}
	if(!Pause){start();	document.getElementById("PausePlay").innerHTML="Pause";}
}

//fillDistanceGrid();
OnChange(0);
//drawCanvas1();
//start();

function SaveData(){
	stop();
	savedData = JSON.stringify({gridSize: gridSize, gridSize2: gridSize2, screenDraw: screenDraw, ColorSpace: ColorSpace, screenDelay: screenDelay, dataGrid: dataGrid});
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

function LoadData(){
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