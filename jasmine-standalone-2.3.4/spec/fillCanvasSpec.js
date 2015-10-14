describe("Canvas Games", function(){

	it("screenDraw should be at zero", function(){
		expect(screenDraw).toEqual(0);
	
	});
	
	it("gridSize should be 4:", function(){
		expect(gridSize).toEqual(4);
	});
	
	it("fixedSpots array should have 4 values", function(){
		setFixedValues();
		
		expect(fixedSpots.length).toEqual(gridSize);
	});
	

});


describe("fixedSpots with gridSize=10",function(){
	beforeEach(function(){
		gridSize = 10;
		setFixedValues();
	});
		
	it("should be points inside the 10x10 grid", function(){
		
		expect(gridSize).toEqual(10);
		expect(fixedSpots[0][0]).toEqual(2);
		expect(fixedSpots[0][1]).toEqual(2);
		expect(fixedSpots[0][2]).toEqual(22);
		expect(fixedSpots[1][0]).toEqual(2);
		expect(fixedSpots[1][1]).toEqual(7);
		expect(fixedSpots[1][2]).toEqual(27);
		expect(fixedSpots[2][0]).toEqual(7);
		expect(fixedSpots[2][1]).toEqual(2);
		expect(fixedSpots[2][2]).toEqual(72);
		expect(fixedSpots[3][0]).toEqual(7);
		expect(fixedSpots[3][1]).toEqual(7);
		expect(fixedSpots[3][2]).toEqual(77);				
	});
	
	
});

describe("rand1",function(){
	beforeEach(function(){
		ColorSpace = 10;
		ColorScale = 10;
	});
	
	it("Should generate values from 0 to ColorSpace",function(){
		var minX = ColorSpace * ColorScale;
		var maxX = -1;
		for (var i=0; i<1000; i++){
			var x = rand1();
			if (x<minX){
				minX = x;
			};
			if (x>maxX){
				maxX = x;
			};
		}
		expect(minX).toBeLessThan(maxX);
		expect(minX).toBeGreaterThan(-1);
		expect(maxX).toBeLessThan(ColorSpace*ColorScale);
	});
});


describe("randGS",function(){
	beforeEach(function(){
		gridSize = 100;
	});
	
	it("Should generate values from 0 to gridSize-1",function(){
		var minX = gridSize;
		var maxX = -1;
		for (var i=0; i<1000; i++){
			var x = randGS(gridSize);
			if (x<minX){
				minX = x;
			};
			if (x>maxX){
				maxX = x;
			};
		}
		expect(minX).toBeLessThan(maxX);
		expect(minX).toBeGreaterThan(-1);
		expect(maxX).toBeLessThan(gridSize);		
	});
});


describe("randPMFr",function(){
	beforeEach(function(){
		forceRange = 10;
	});
	
	it("Should generate values from -forceRange to +forceRange",function(){
		var minX = forceRange;
		var maxX = -forceRange;
		for (var i=0; i<1000; i++){
			var x = randPMFr();
			if (x<minX){
				minX = x;
			};
			if (x>maxX){
				maxX = x;
			};
		}
		expect(minX).toEqual(-forceRange);
		expect(maxX).toEqual(forceRange);		
	});
});

describe("randPMZ1",function(){
	
	it("Should generate values -1,0,1",function(){
		var negX = false;
		var posX = false;
		var zeroX = false;
		var outOfBounds = false;
		
		for (var i=0; i<1000; i++){
			var x = randPZM1();
			switch(x){
				case 0: {
					zeroX = true;
					break;
				}
				case 1: {
					posX = true;
					break;
				}
				case -1: {
					negX = true;
					break;
				}
				default:{
					outOfBounds = true;
					break;
				}
			}
		};
		expect(posX).toEqual(true);
		expect(negX).toEqual(true);		
		expect(zeroX).toEqual(true);
		expect(outOfBounds).toEqual(false);

	});
});

describe("grid cell initialization",function(){
	var x =  new gridCell(5,4,3);
	
	it("gridcell should be 5,4,3",function(){
		expect(x.r).toEqual(5);
		expect(x.g).toEqual(4);
		expect(x.b).toEqual(3);
	});
	
});

describe("fillDistanceGrid",function(){
	beforeEach(function(){
		gridSize = 10;
		gridSize2 = gridSize*gridSize;
		distanceGrid = Array.matrix(gridSize,gridSize,0);
		fillDistanceGrid();
	});
	
	it("check that the distance grid is filled out",function(){
		expect(distanceGrid[0][0].dist).toEqual(1000);
		expect(distanceGrid[0][0].dist2).toEqual(1000000);
		expect(distanceGrid[gridSize-1][0].dist).toEqual(9);
		expect(distanceGrid[0][gridSize-1].dist).toEqual(9);
		expect(distanceGrid[gridSize-1][0].dist2).toEqual(81);
		expect(distanceGrid[0][gridSize-1].dist2).toEqual(81);
		expect(Math.floor(distanceGrid[gridSize-1][gridSize-1].dist)).toEqual(12);
		expect(distanceGrid[gridSize-1][gridSize-1].dist2).toEqual(162);
		expect(Math.floor(distanceGrid[5][5].dist)).toEqual(7);
		expect(distanceGrid[5][5].dist2).toEqual(50);
	});
	
});
