// C Code for FFT: http://fourier.eng.hmc.edu/e101/lectures/Image_Processing/node5.html

// RGB data structure
function gridCell(rc,gc,bc){this.r=rc;this.g=gc;this.b=bc;}

// element swap function
// function SWAP(a,b){
// 	var t = a;
// 	a = b;
// 	b = t;
// }

// 1-D FFT for a line of data
function fftLine(xr, xi, N){

  var m = Math.ceil(Math.log2(N));
  for (var i=0; i<N; ++i) {        // bit reversal 
    var j=0;
    for (var k=0; k<m; ++k)
      j=(j << 1) | (1 & (i >> k));
	    if (j < i) 
    	  { var t = xr[i]; xr[j]=xr[i]; xr[i]=t; t=xi[i]; xi[i]=xi[j]; xi[j]=t; } // swap
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
			var tempr=xr[j1+n]*c-xi[j1]*s;
			var tempi=xi[j1+n]*c+xr[j1]*s;
			xr[j1+n]=xr[j1]-tempr;
			xi[j1+n]=xi[j1]-tempi;
			xr[j1]=xr[j1]+tempr;
			xi[j1]=xi[j1]+tempi;
		}
      }
      k+=2*n;
    }
  }
  
}



// function call to calculate 2D FFT
function fft2d(img2dRGB, N){
    var M = Math.pow(2,Math.ceil(Math.log2(N)));
	var N2 = N*N;
	var v = 256; // maximum video range
	var fft = new Array(N2);
	
	for (var i = 0; i<N2; i++){
		fft[i] = new gridCell(0,0,0);
	}
	var lineR = new Array(N);
	var lineRi = new Array(N);
	var lineG = new Array(N);
	var lineGi = new Array(N);
	var lineB = new Array(N);
	var lineBi = new Array(N);
	// fft by row and color component
	for (var row = 0; row< N; row++){
		var kR = row * N;
		// copy out line of data...
		for (var col = kR; col < kR + M; col++){
			var kc = col-kR; // convenience column number
			if(kc<N){
				lineR[kc] = img2dRGB[col].r;
				lineRi[kc] = 0;
				lineG[kc] = img2dRGB[col].g;
				lineGi[kc] = 0;
				lineB[kc] = img2dRGB[col].b;
				lineBi[kc] = 0;
			} else {
				lineR[kc] = 0;
				lineRi[kc] = 0;
				lineG[kc] = 0;
				lineGi[kc] = 0;
				lineB[kc] = 0;
				lineBi[kc] = 0;
			}
		}
		// transform each line
		fftLine(lineR, lineRi, M);
		fftLine(lineG, lineGi, M);
		fftLine(lineB, lineBi, M);
		var lineRMag = [];
		var lineGMag = [];
		var lineBMag = [];
		// normalize each line Z
		for ( var i =0; i<M; i++){
			lineRMag[i] = Math.sqrt(lineR[i]*lineR[i]+lineRi[i]*lineRi[i]);
			lineGMag[i] = Math.sqrt(lineG[i]*lineG[i]+lineGi[i]*lineGi[i]);
			lineBMag[i] = Math.sqrt(lineB[i]*lineB[i]+lineBi[i]*lineBi[i]);
		}
		var lineRMax = v/(Math.max(...lineRMag) + 0.0000001);
		var lineGMax = v/(Math.max(...lineGMag) + 0.0000001);
		var lineBMax = v/(Math.max(...lineBMag) + 0.0000001);
		// write data out to fft array
		for (var col = kR; col < kR + N; col++){
			var kc2 = Math.floor(col-kR-(N-1)/2+M)%M; // convenience variable splitting spectrum in half 
			fft[col].r = Math.floor(lineRMag[kc2]*lineRMax);
			fft[col].g = Math.floor(lineGMag[kc2]*lineGMax);
			fft[col].b = Math.floor(lineBMag[kc2]*lineBMax);
		}
		
	}




return fft;



}