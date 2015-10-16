// C Code for FFT: http://fourier.eng.hmc.edu/e101/lectures/Image_Processing/node5.html

// RGB data structure
function gridCellC(rc,rci,gc,gci,bc,bci){this.r=rc;this.ri=rci;this.g=gc;this.gi=gci;this.b=bc;this.bi=bci;}

// element swap function
// function SWAP(a,b){
// 	var t = a;
// 	a = b;
// 	b = t;
// }

// 1-D FFT for a line of data
// function fftLine(xr, xi, N){
// 
//   var m = Math.ceil(Math.log2(N));
//   for (var i=0; i<N; ++i) {        // bit reversal 
//     var j=0;
//     for (var k=0; k<m; ++k)
//       j=(j << 1) | (1 & (i >> k));
// 	    if (j < i) 
//     	  { var t = xr[i]; xr[j]=xr[i]; xr[i]=t; t=xi[i]; xi[i]=xi[j]; xi[j]=t; } // swap
//   }
//   for (var i=0; i<m; i++) {         // for log N stages 
//     var n=Math.pow(2.0,i);  
//     var w=Math.PI/n;
//     //if (inverse) w=-w;
//     var k=0;
//     while (k<N-1) {             // for N components 
//       for (var j=0; j<n; j++) {     // for each section 
// 		var c=Math.cos(-j*w); var s=Math.sin(-j*w); 
// 		var j1=k+j;
// 		if(j1+n<N){
// 			var tempr=xr[j1+n]*c-xi[j1]*s;
// 			var tempi=xi[j1+n]*c+xr[j1]*s;
// 			xr[j1+n]=xr[j1]-tempr;
// 			xi[j1+n]=xi[j1]-tempi;
// 			xr[j1]=xr[j1]+tempr;
// 			xi[j1]=xi[j1]+tempi;
// 		}
//       }
//       k+=2*n;
//     }
//   }
// }

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


// fold and clip the FFT array to the middle NxN region
function foldAndClipArray(dataC, N, M){
    var M2 = M*M;
	var N2 = N*N;
	var midM = M/2;
	var midN = N/2;
	
	// initialize empty complex data array for R+G+B
	var fCA = new Array(N2);
	for (var i = 0; i<N2; i++){
		var kRN = Math.floor(i/N);
		var kRM = (kRN+midM+midN)%M;
		var kCN = (i-kRN*N);
		var kM = kRM*M+(kCN+midM+midN)%M;
		fCA[i] = new gridCellC(dataC[kM].r,0,dataC[kM].g,0,dataC[kM].b,0);
	}
	return fCA;
}

// function call to calculate 2D FFT
function fft2d(img2dRGB, N){
    var M = Math.pow(2,Math.ceil(Math.log2(N)));
    var M2 = M*M;
	var N2 = N*N;
	var v = 256; // maximum video range
	var fft = moveRealToComplex(img2dRGB,N,M);
		
	// var lineR = new Array(N);
// 	var lineRi = new Array(N);
// 	var lineG = new Array(N);
// 	var lineGi = new Array(N);
// 	var lineB = new Array(N);
// 	var lineBi = new Array(N);
	// fft by row and color component - pull out multiple lines of data...
	for (var row = 0; row< M; row++){
	// 	var kR = row * N;
// 		// copy out line of data...
// 		for (var col = kR; col < kR + M; col++){
// 			var kc = col-kR; // convenience column number
// 			if(kc<N && kR<N){
// 				lineR[kc] = img2dRGB[col].r;
// 				lineRi[kc] = 0;
// 				lineG[kc] = img2dRGB[col].g;
// 				lineGi[kc] = 0;
// 				lineB[kc] = img2dRGB[col].b;
// 				lineBi[kc] = 0;
// 			} else {
// 				lineR[kc] = 0;
// 				lineRi[kc] = 0;
// 				lineG[kc] = 0;
// 				lineGi[kc] = 0;
// 				lineB[kc] = 0;
// 				lineBi[kc] = 0;
// 			}
// 		}
		
		
		// transform each line
	// 	fftLine(lineR, lineRi, M);
// 		fftLine(lineG, lineGi, M);
// 		fftLine(lineB, lineBi, M);
		
		var line = extractRow(fft,M,row,"r");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"r");
		line = extractRow(fft,M,row,"g");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"g");
		line = extractRow(fft,M,row,"b");
		fftLineC(line,M);
		insertRow(line,fft,M,row,"b");
		
		// copy FFT back to temporary array
	// 	var kRM = row*M; // convenience row pointer
// 		for (var col=kRM; col<kRM+M; col++){
// 			var kCM = col-kRM;
// 			fft[col].rc = lineR[kCM];
// 			fft[col].rci = lineRi[kCM];
// 			fft[col].gc = lineG[kCM];
// 			fft[col].gci = lineGi[kCM];
// 			fft[col].bc = lineB[kCM];
// 			fft[col].bci = lineBi[kCM];
// 		}		
	}


	var fftReturn = foldAndClipArray(fft, N, M);

	return fftReturn;



}









// 
// 		var lineRMag = [];
// 		var lineGMag = [];
// 		var lineBMag = [];
// 		// normalize each line Z
// 		for ( var i =0; i<M; i++){
// 			lineRMag[i] = Math.sqrt(lineR[i]*lineR[i]+lineRi[i]*lineRi[i]);
// 			lineGMag[i] = Math.sqrt(lineG[i]*lineG[i]+lineGi[i]*lineGi[i]);
// 			lineBMag[i] = Math.sqrt(lineB[i]*lineB[i]+lineBi[i]*lineBi[i]);
// 		}
// 		var lineRMax = v/(Math.max(...lineRMag) + 0.0000001);
// 		var lineGMax = v/(Math.max(...lineGMag) + 0.0000001);
// 		var lineBMax = v/(Math.max(...lineBMag) + 0.0000001);
// 		// write data out to fft array
// 		for (var col = kR; col < kR + N; col++){
// 			var kc2 = Math.floor(col-kR-(N-1)/2+M)%M; // convenience variable splitting spectrum in half 
// 			fft[col].r = Math.floor(lineRMag[kc2]*lineRMax);
// 			fft[col].g = Math.floor(lineGMag[kc2]*lineGMax);
// 			fft[col].b = Math.floor(lineBMag[kc2]*lineBMax);
// 		}