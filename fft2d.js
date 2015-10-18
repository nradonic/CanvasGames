// C Code for FFT: http://fourier.eng.hmc.edu/e101/lectures/Image_Processing/node5.html

// RGB data structure
function gridCellC(rc,rci,gc,gci,bc,bci){this.r=rc;this.ri=rci;this.g=gc;this.gi=gci;this.b=bc;this.bi=bci;}


// 1-D FFT for a line of data
function fftLineC(x, N){

  var m = Math.ceil(Math.log2(N));
  for (var i=0; i<N; ++i) {        // bit reversal 
    var j=0;
    for (var k=0; k<m; ++k)
      j=(j << 1) | (1 & (i >> k));
	    if (j < i) 
    	  { var t = x[i].r; x[j].r=x[i].r; x[i].r=t; t=x[i].i; x[i].i=x[j].i; x[j].i=t; } // swap
  }
  for (var i=0; i<m; i++) {         // for log N stages 
    var n=Math.pow(2.0,i);  
    var w=Math.PI/n;
    //if (inverse) w=-w;
    var k=0;
    while (k<N-1) {             // for N components 
      for (var j=0; j<n; j++) {     // for each section 
		var c=Math.cos(-j*w); var s=Math.sin(-j*w); 
		var j1=k+j;
		if(j1+n<N){
			var tempr=x[j1+n].r*c-x[j1].i*s;
			var tempi=x[j1+n].i*c+x[j1].r*s;
			x[j1+n].r=x[j1].r-tempr;
			x[j1+n].i=x[j1].i-tempi;
			x[j1].r=x[j1].r+tempr;
			x[j1].i=x[j1].i+tempi;
		}
      }
      k+=2*n;
    }
  }
}

// transfer data from smaller REAL RGB array to larger 2^M complex array
function moveRealToComplex(dataR, N, M){
    var M2 = M*M;
	var N2 = N*N;
	
	// initialize empty complex data array for R+G+B
	var fft = new Array(M2);
	for (var i = 0; i<M2; i++){
		fft[i] = new gridCellC(0,0,0,0,0,0);
	}
	//copy real data into complex array
	for(var row=0; row<M; row++){
		var kRM=row*M; // convenience variable
		var kRN=row*N; // convenience variable
		// copy out line of data...
		for (var col = kRM; col < kRM + M; col++){
			var kCol=col-kRM;
			var kCN = kRN + kCol; // convenience column number
			if(kCol<N && row<N){
				fft[col].r = dataR[kCN].r;
				fft[col].g = dataR[kCN].g;
				fft[col].b = dataR[kCN].b;
			} else {
				fft[col].r = 0;
				fft[col].g = 0;
				fft[col].b = 0;
			}
			fft[col].ri = 0;
			fft[col].gi = 0;
			fft[col].bi = 0;
		}
	}
	return fft;
}

// pull out the needed row of data for that color
function extractRow(fftData,M,row, color){
	var line = new Array(M);
	var kR = row * M;
	// copy out line of data...
	for (var col = kR; col < kR + M; col++){
		var kc = col-kR; // convenience column number
		switch(color){
		case "r":{
					line[kc]={
								r: fftData[col].r,
								i: fftData[col].ri
							}
					break;
				}
		case "g":{
					line[kc]={
								r: fftData[col].g,
								i: fftData[col].gi
							}
					break;
				}
		case "b":{
					line[kc]={
								r: fftData[col].b,
								i: fftData[col].bi
							}
					break;
				}
		}
	}
	return line;
}

// pull out the needed row of data for that color
function extractColumn(fftData,M,column, color){
	var M2 = M*M;
	var line = new Array(M);
	// copy out line of data...
for (var cell = 0; cell<M; cell++){
		var kCM = cell*M+column;
		switch(color){
		case "r":{
					line[cell]={
								r: fftData[kCM].r,
								i: fftData[kCM].ri
							}
					break;
				}
		case "g":{
					line[cell]={
								r: fftData[kCM].g,
								i: fftData[kCM].gi
							}
					break;
				}
		case "b":{
					line[cell]={
								r: fftData[kCM].b,
								i: fftData[kCM].bi							}
					break;
				}
		}
	}
	return line;
}

// put the complex data line back into the fft array
function insertRow(line,fftData,M,row,color){
	var kM = row*M;
	for (var iX = kM; iX<kM+M; iX++){
		var kCM = iX - kM;
		switch(color){
		case "r":{
					fftData[iX].r = line[kCM].r;
					fftData[iX].ri = line[kCM].i;
					break;
				}
		case "g":{
					fftData[iX].g = line[kCM].r;
					fftData[iX].gi = line[kCM].i;
					break;
				}
		case "b":{
					fftData[iX].b = line[kCM].r;
					fftData[iX].bi = line[kCM].i;
					break;
				}
		}
	}
}

// put the complex data line back into the fft array
function insertColumn(line,fftData,M,column,color){
	var M2 = M*M;
	for (var iX = 0; iX<M; iX++){
		var kCM = column + iX*M;
		switch(color){
		case "r":{
					fftData[kCM].r = line[iX].r;
					fftData[kCM].ri = line[iX].i;
					break;
				}
		case "g":{
					fftData[kCM].g = line[iX].r;
					fftData[kCM].gi = line[iX].i;
					break;
				}
		case "b":{
					fftData[kCM].b = line[iX].r;
					fftData[kCM].bi = line[iX].i;
					break;
				}
		}
	}
}

// fold and clip the FFT array to the middle NxN region
function foldAndClipArray(dataC, N, M){
    var M2 = M*M;
	var N2 = N*N;
	var midM = Math.floor(M/2);
	var midN = Math.floor(N/2);
	
	// initialize empty complex data array for R+G+B
	var fCA = new Array(N2);
	for (var i = 0; i<N2; i++){
		var kRN = Math.floor(i/N);
		var kRM = (kRN+M-midN)%M;
		var kCN = (i-kRN*N);
		var kM = kRM*M+(kCN+M-midN)%M;
		fCA[i] = new gridCellC(dataC[kM].r,0,dataC[kM].g,0,dataC[kM].b,0);
		if(kM===0){ fCA[i].r=0; fCA[i].g=0;fCA[i].b=0;}
	}
	return fCA;
}


// normalize fft results to 256
function normalize(fftData, N){
	var N2 = N*N;
	var mr = 0;
	var midN = Math.floor(N/2);
	var midCell = N*midN+midN;
	for (var i=0; i<N2; i++){
		var r = fftData[i].r;
		var ri = fftData[i].ri;
		var g = fftData[i].g;
		var gi = fftData[i].gi;
		var b = fftData[i].b;
		var bi = fftData[i].bi;
		if(i!==midCell){
			mr = Math.max(mr,r*r+ri*ri,g*g+gi*gi,b*b+bi*bi);
		}
	}
	mr = 256/(Math.sqrt(mr+.0000001));
	for (var i=0; i<N2; i++){
		var r = fftData[i].r;
		var ri = fftData[i].ri;
		var g = fftData[i].g;
		var gi = fftData[i].gi;
		var b = fftData[i].b;
		var bi = fftData[i].bi;
		fftData[i].r = Math.floor(Math.sqrt(r*r+ri*ri)*mr);
		fftData[i].g = Math.floor(Math.sqrt(g*g+gi*gi)*mr);
		fftData[i].b = Math.floor(Math.sqrt(b*b+bi*bi)*mr);
	}
	fftData[midCell].r = 255;
	fftData[midCell].g = 255;
	fftData[midCell].b = 255;
}


// function call to calculate 2D FFT
function fft2d(img2dRGB, N){
    var M = Math.pow(2,Math.ceil(Math.log2(N)));
    var M2 = M*M;
	var N2 = N*N;
	var v = 256; // maximum video range
	var fft = moveRealToComplex(img2dRGB,N,M);
		
	// fft by row and color component - pull out multiple lines of data...
	for (var row = 0; row< M; row++){
		
		// transform each line		
		var line = extractRow(fft,M,row,"r");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"r");
		
		line = extractRow(fft,M,row,"g");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"g");
		
		line = extractRow(fft,M,row,"b");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"b");	
	}
	
	for (var col=0; col<M; col++){
		//transform each column
		var line = extractColumn(fft,M,col,"r");
		fftLineC(line,M);
		insertColumn(line,fft,M,col,"r");
		
		line = extractColumn(fft,M,col,"g");
		fftLineC(line,M);
		insertColumn(line,fft,M,col,"g");
		
		line = extractColumn(fft,M,col,"b");
		fftLineC(line,M);
		insertColumn(line,fft,M,col,"b");	
	}

	var fftReturn = foldAndClipArray(fft, N, M);
	normalize(fftReturn, N);
	
	return fftReturn;



}





