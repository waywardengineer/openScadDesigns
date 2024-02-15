sides = 10;
angleIncrementPercentage = 100;
radius = 25;
lineWidth = 0.48;
shapeType = 1;
branchAtCCWsegments = [7, 15, 0, 0];
branchAtCWsegments = [7, 15, 0, 0];
branchfterLevelsCCW = [0, 0, 0, 0];
branchfterLevelsCW = [0, 0, 0, 0];
stopAtLastBranch = true;
polygonOffsetPercentage = 0;
polygonOffsetAmountX100 = 0;
shapeOffsetAmountx100 = 0;
polygonOffsetAmount = polygonOffsetAmountX100 / 100;
shapeOffsetAmount = shapeOffsetAmountx100 / 100;

startAtBranchLevel = 0;
stopAtBranchLevel = 0;
rotationalCopies = 1;

stopAtLevels = [stopAtLastBranch ? max(branchAtCCWsegments) : 0, stopAtLastBranch ? max(branchAtCWsegments) : 0];
branchAtSegments = [branchAtCCWsegments, branchAtCWsegments];
branchAfterLevels = [branchfterLevelsCCW, branchfterLevelsCW];
angleIncrement = (180/sides) * (angleIncrementPercentage / 100);

module regular_polygon(sides, radius, wallThickness = 0){
     angles=[ for (i = [0:sides-1]) i*(360/sides) ];
     coords=[ for (th=angles) [radius*cos(th), radius*sin(th)] ];
     if (wallThickness <= 0){
         polygon(coords);
     }
     else {
        difference(){
            polygon(coords);
            offset(r=-wallThickness) polygon(coords);
        }
    }
}
module triangle(
    points, 
    direction = 1, 
    lineWidth = lineWidth, 
    shapeType = shapeType, 
    polygonOffsetPercentage = polygonOffsetPercentage, 
    polygonOffsetAmount = polygonOffsetAmount
    ){
    polygonOffset =  abs(polygonOffsetAmount) > 0 ? polygonOffsetAmount : polygonOffsetPercentage * maxTriangleSide(points[0], points [1], points[2]) / 100; 
    if (shapeType == 0){
        offsetLineTriangle (points, lineWidth, -polygonOffsetAmount, direction);
    } else if (shapeType == 1){
        offset(r=polygonOffset) polygon(points);
    } else if (shapeType == 2){     
        notchedOffsetTriangle(points, lineWidth * 1.5);
    } else if (shapeType == 3){     
        notchedOffsetTriangle(points, lineWidth * 0.5);
    } else if (shapeType == 4){
        circleInTriangle(points,  circleRatio = 3);
    } else if (shapeType == 5){
        difference(){
            offset(r=-0.1) polygon(points);
            circleInTriangle(points,  circleRatio = 1);
        }
    } 
    
}

module notchedOffsetTriangle(points, outerOffset=0, circleSize = 1, pointToNotch = 1){
    p1 = points[0];
    p2 = points[1];
    p3 = points[2];
    lengths = triangleLengths(p1, p2, p3);
    angleC = angleFromTriangleSides(lengths[2], lengths[0], lengths[1]);
    angleB = angleFromTriangleSides(lengths[1], lengths[0], lengths[2]);
    a1 = lengths[0] / (tan(angleC/2)/tan(angleB/2) + 1);
    maxOffset = tan(angleC/2) * a1;
    if (outerOffset < maxOffset) {
        a1Angle = angleOfVector(xdiff(p2, p3), ydiff(p2, p3));
        baseIntersectionPoint = subtractPoints(p3, polarToCartesian(a1, 90 -a1Angle));
        wonkyVector = polarToCartesian(maxOffset,   0- a1Angle );
        stupidVector = [
            p1[0] <  baseIntersectionPoint[0] ? -abs(wonkyVector[0]): abs(wonkyVector[0]),
            p1[1] <  baseIntersectionPoint[1] ? -abs(wonkyVector[1]): abs(wonkyVector[1]), 
        ];
        center = addPoints (baseIntersectionPoint, stupidVector );
        offsetRatio = 1 - outerOffset / maxOffset;
        newPoints = [
            for (i = [0:2])
                addPoints(center, [offsetRatio * xdiff(center, points[i]), offsetRatio * ydiff(center, points[i])])
            
        ];
        difference(){
            polygon(newPoints);
            translate(newPoints[pointToNotch]) circle(circleSize);
        }
    }
}

module circleInTriangle(points, circleRatio = 1, locationType = 0){
    p1 = points[0];
    p2 = points[1];
    p3 = points[2];
    lengths = triangleLengths(p1, p2, p3);
    angleA = angleFromTriangleSides(lengths[0], lengths[2], lengths[1]);
    angleC = angleFromTriangleSides(lengths[2], lengths[0], lengths[1]);
    angleB = angleFromTriangleSides(lengths[1], lengths[0], lengths[2]);
    a1 = lengths[0] / (tan(angleC/2)/tan(angleB/2) + 1);
    distanceToCenter = tan(angleC/2) * a1;
    a1Angle = angleOfVector(xdiff(p2, p3), ydiff(p2, p3));
    aMidPoint = subtractPoints(p3, polarToCartesian(a1, 90 -a1Angle));
    wonkyVector = polarToCartesian(distanceToCenter,   0- a1Angle );
    stupidVector = [
        p1[0] <  aMidPoint[0] ? -abs(wonkyVector[0]): abs(wonkyVector[0]),
        p1[1] <  aMidPoint[1] ? -abs(wonkyVector[1]): abs(wonkyVector[1]), 
    ];
    center = addPoints (aMidPoint, stupidVector );
    p1p3Mid = addPoints (p1, [0.5 * xdiff(p1, p3), 0.5 * ydiff(p1, p3)]);
    if (locationType == 0){
        translate(center) circle(distanceToCenter * circleRatio);       
    } else if (locationType == 1){
        steps = 5;
        for (i = [0:steps]){
            translate(interpolatePoints(p1, p3, i / steps)) circle(distanceToCenter * circleRatio);
        }
        
    }

}


module recursivePolygonSide(
        sides = sides, 
        rotation = 0, 
        angleIncrement = angleIncrement, 
        radius = radius, 
        minRadius = 2, 
        segmentLevel = 0,
        branchLevel = 0,
        startAtBranchLevel = startAtBranchLevel,
        stopAtBranchLevel = stopAtBranchLevel,
        shapeType = shapeType, 
        branchAtSegments = branchAtSegments,
        stopAtLevels = stopAtLevels,
        polygonOffsetPercentage = polygonOffsetPercentage, 
        polygonOffsetAmount = polygonOffsetAmount,
    ) {
    innerRadius = radius * polygonSizeRatio(sides, abs(angleIncrement));
    direction = angleIncrement < 0 ? -1: 1;
    p2PolygonRotation = direction * (360 / sides);
    p1 = polarToCartesian(radius, rotation);
    p2 = polarToCartesian(innerRadius, rotation + angleIncrement - p2PolygonRotation); 
    p3 = polarToCartesian(innerRadius, rotation + angleIncrement);
    branchingParamIndex = angleIncrement < 0 ? 0: 1;
    newRadius = innerRadius;
    if (stopAtLevels[branchingParamIndex] == 0 || segmentLevel <= stopAtLevels[branchingParamIndex]){
        if (startAtBranchLevel <= branchLevel && (stopAtBranchLevel == 0 || branchLevel < stopAtBranchLevel)){
            triangle(
                [p1, p2, p3], 
                direction = direction, 
                shapeType = shapeType, 
                polygonOffsetPercentage = polygonOffsetPercentage, 
                polygonOffsetAmount = polygonOffsetAmount
            );
        }

        for (i = [0:3]){
            branchAtSegment = branchAtSegments[branchingParamIndex][i];
            branchAfterLevel = branchAfterLevels[branchingParamIndex][i];
            if (branchAtSegment > 0 && segmentLevel == branchAtSegment && branchLevel >= branchAfterLevel){
                lengths = triangleLengths(p1, p2, p3);
                branchRadius = radius * lengths[1] / lengths[2];
                branchTranslationToOrigin = polarToCartesian(branchRadius, rotation);
                rotate (180)
                translate ([-branchTranslationToOrigin[0] -p3[0], -branchTranslationToOrigin[1] - p3[1]])
                    recursivePolygonSide(
                        rotation = rotation, 
                        angleIncrement = -angleIncrement, 
                        radius = branchRadius, 
                        segmentLevel = 0,
                        branchLevel = branchLevel + 1, 
                        shapeType = shapeType, 
                        branchAtSegments = branchAtSegments,
                        polygonOffsetPercentage = polygonOffsetPercentage, 
                        polygonOffsetAmount = polygonOffsetAmount
                    );
            }
        }
        if(newRadius > minRadius){
            recursivePolygonSide(
                rotation = rotation + angleIncrement, 
                angleIncrement = angleIncrement, 
                radius = newRadius, 
                segmentLevel = segmentLevel + 1, 
                branchLevel = branchLevel, 
                shapeType = shapeType, 
                branchAtSegments = branchAtSegments,
                polygonOffsetPercentage = polygonOffsetPercentage, 
                polygonOffsetAmount = polygonOffsetAmount
            );
        }
    }
}

module offsetLineTriangle (points, lineSpacing, outerOffset, direction = 1){
    p1 = points[0];
    p2 = points[1];
    p3 = points[2];
    baseAngle = angleOfVector(xdiff(p1, p2), ydiff(p1, p2));
    lengths = triangleLengths (p1, p2, p3);
    topInternalAngle = angleFromTriangleSides(lengths[0], lengths[1], lengths[2]);
    extensionAngle = 180 - topInternalAngle;
    extensionLength = abs(lengths[2] * cos(extensionAngle));
    totallineLength = lengths[0] + extensionLength;
    lineAreaSpan = abs(lengths[2] * sin(extensionAngle));
    numLineSubtractions = lineAreaSpan / (2 * lineSpacing);
    intersection(){
        offset(r=-outerOffset) polygon(points);
        difference (){
            polygon(points);
            for (i = [0:numLineSubtractions]){
                translate (p1) 
                    rotate (-baseAngle)  
                    translate([direction * ((1 + i * 2 + (direction < 0 ? 1 : 0)) * lineSpacing + outerOffset ), -extensionLength]) 
                    square([lineSpacing, totallineLength], false);
            }
        }
    }
}



function xdiff(p1, p2) = p2[0] - p1[0];

function ydiff(p1, p2) = p2[1] - p1[1];

function lineLength (p1, p2) = sqrt(xdiff(p1, p2) ^ 2 + ydiff(p1, p2)^2);

function maxTriangleSide (pa, pb, pc) = max(lineLength(pa, pb), lineLength(pb, pc), lineLength(pc, pa));

function angleFromTriangleSides(op, aj1, aj2) = acos ((aj1 ^ 2 + aj2 ^ 2 - op ^ 2) / ( 2 * aj1 * aj2));

function edgeCaseAtanAdjustment(x, y) = x == 0 ? (y < 0 ? 180 : 0) : (x < 0 ? 180 : 0);

function polarToCartesian(r, theta) = [r * cos(theta), r * sin(theta)];

function angleOfVector (x, y) = (x == 0 ? 0 : 90 - atan(y / x)) + edgeCaseAtanAdjustment(x, y);

function triangleLengths(pa, pb, pc) = [lineLength(pb, pc), lineLength(pa, pc), lineLength(pa, pb)];

function sideIntersectionRatio(polygonSides = sides, angleIncrement = angleIncrement) = sin(polygonAngle - angleIncrement) / sin(angleIncrement);

function polygonSizeRatio(polygonSides = sides, angleIncrement = angleIncrement) = 
    sin(180 - (360 / polygonSides)) / ((1 + sideIntersectionRatio(polygonSides, angleIncrement)) *  sin(angleIncrement));

function addPoints (p1, p2) = [p1[0] + p2[0], p1[1] + p2[1]];

function subtractPoints (p1, p2) = [p1[0] - p2[0], p1[1] - p2[1]];

function interpolatePoints (p1, p2, ratio) = addPoints (p1, [ratio * xdiff(p1, p2), ratio * ydiff(p1, p2)]);
polygonAngle = 360 / sides;
sizeRatio = polygonSizeRatio(sides, angleIncrement);
angles=[ for (i = [0:sides-1]) i*(360/sides)];
everyOtherAngle = [for (i = [0:3:sides-1]) angles[i]];    
halfPolygonAngle = 180 / sides;
outerPolygonInnerAngles = [for (angle = angles) angle + 180 + halfPolygonAngle];
verticalOffset = radius * sin(halfPolygonAngle) * 2 / (1 + sideIntersectionRatio());
//recursivePolygonSide ();

copiesToMake = rotationalCopies - 1;
copyAngleDivision = 360 / (rotationalCopies);

for (i=[0:copiesToMake]){
    rotate(i * copyAngleDivision){
        translate([radius * 2 * cos(360/(2*sides)), verticalOffset])
        rotate(outerPolygonInnerAngles[0])
            //difference(){
                offset(r=shapeOffsetAmount)
                recursivePolygonSide ();
                //recursivePolygonSide (polygonOffsetAmount = -0.5);
            //}
        rotate(halfPolygonAngle) 
            offset(r=shapeOffsetAmount)
            recursivePolygonSide (branchAtSegments = [[0, 0, 0, 0], [0, 0, 0, 0]]);
    }
}



