sides = 6;
angleIncrement = 1;
radius = 20;
minRadius = 3;
heightPerLayer = 2;
twist = 0;
inverse = true;
maxRadius = 50;
fixedRadius = 30;
sinRadius = 20;
height = 100;
cornerRadius = 5;
polygonAngle = 360 / sides;
sideIntersectionRatio =  sin(polygonAngle - angleIncrement) / sin(angleIncrement);
sizeRatio = sin(180 - polygonAngle) / ((1 + sideIntersectionRatio) * sin(angleIncrement));
paths = [[ for (i = [0:sides-1]) i], [ for (i = [0:sides-1]) i + sides]];
   


module regular_polygon(sides, radius){
     angles=[ for (i = [0:sides-1]) i*(360/sides) ];
     coords=[ for (th=angles) [radius*cos(th), radius*sin(th)] ];
     polygon(coords);
}
 
module recursivePolygon(sides, angles, angleIncrement, radius, sizeRatio, minRadius, paths) {
    innerAngles = [for (a = angles) each a + angleIncrement];
    innerRadius = radius * sizeRatio;
    outerCoords = [for (th=angles) [radius*cos(th), radius*sin(th)]];
    innerCoords = [for (th=innerAngles) [innerRadius*cos(th), innerRadius*sin(th)]];
    polygon([each outerCoords, each innerCoords], paths);
    newRadius = innerRadius * sizeRatio;
    if(newRadius > minRadius){
        recursivePolygon(sides, [for (a = innerAngles) each a + angleIncrement], angleIncrement, newRadius, sizeRatio, minRadius, paths);
    }
}
module recursivePolygonSinSegment(angles, radius, currentHeight) {
    innerAngles = [for (a = angles) each a + angleIncrement];
    //innerRadius = radius * sizeRatio;
   
    outerCoords = [for (th=angles) [radius*cos(th), radius*sin(th)]];
    //innerCoords = [for (th=innerAngles) [innerRadius*cos(th), innerRadius*sin(th)]];
    linear_extrude(heightPerLayer, twist = twist)
        polygon(points = outerCoords);
    newHeight = currentHeight + heightPerLayer;
    if(newHeight <= height){
        newRadius = fixedRadius + sinRadius * sin(180*(newHeight/height));
        translate([0, 0, heightPerLayer])
            recursivePolygonSinSegment([for (a = innerAngles) each a + angleIncrement], newRadius, newHeight);
    }
    if(newHeight > height){
        newRadius = radius * sizeRatio ^ 2 ;
        translate([0, 0, heightPerLayer])
            recursivePolygonStraightSegment([for (a = innerAngles) each a + angleIncrement], newRadius);
    }}

module recursivePolygonStraightSegment(angles, radius) {
    innerAngles = [for (a = angles) each a + angleIncrement];
    outerCoords = [for (th=angles) [radius*cos(th), radius*sin(th)]];
    linear_extrude(heightPerLayer, twist = twist)
        polygon(points = outerCoords);
    newRadius = radius * sizeRatio ^ 2 ;
    if(newRadius > minRadius && newRadius < maxRadius){
        translate([0, 0, heightPerLayer])
            recursivePolygonStraightSegment([for (a = innerAngles) each a + angleIncrement], newRadius);
    }
}

module layer(points, height) {
    linear_extrude(heightPerLayer, twist = twist)
        translate([-points[0][0], 0])
        minkowski() {
            polygon(points = points);
            translate(points[0])
            circle(cornerRadius);
        }
}


module recursivePolygonBase(){
    angles=[ for (i = [0:sides-1]) i*(360/sides) ];
    outerCoords = [for (th=angles) [radius*cos(th), radius*sin(th)]];
    recursivePolygonSinSegment (angles, fixedRadius, 0); 
        
}


recursivePolygonBase();