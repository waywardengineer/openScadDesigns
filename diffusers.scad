
sides = 8;
layerHeight = 1;
fixedRadius = 30;
sinRadius = 15;
twistPer100mm = 150;
segmentEndPoints = [60, 160];
segmentPatternHeights = [60, 100];
hullInterval = 2;
twist = 0;

cornerRadius = 0.6;
cornerLength = 2;
cornerVariableInset =  0;
cornerFixedInset = 1;//cornerRadius - 1;
cornerFixedOffsetAngle = 0;
cornerVariableOffsetAngle = -2;

segmentHeights = [for (i = [0:len(segmentEndPoints)-1]) segmentEndPoints[i] - (i == 0 ? 0 : segmentEndPoints[i-1])]; 
segmentSteps = [for (segmentHeight=segmentHeights) segmentHeight / layerHeight];
segmentEndSteps = [for (segmentEndPoint=segmentEndPoints) segmentEndPoint / layerHeight];
segmentPatternSteps = [for (segmentPatternHeight=segmentPatternHeights) segmentPatternHeight / layerHeight];

polygonAngle = 360 / sides;
angleIncrement = twistPer100mm * (layerHeight / 100);
sideIntersectionRatio =  sin(polygonAngle - angleIncrement) / sin(angleIncrement);
sizeRatio = sin(180 - polygonAngle) / ((1 + sideIntersectionRatio) * sin(angleIncrement));
angles=[for (j = [0:sides-1]) j*(360/sides) ];


module corner(radius, height) {
    linear_extrude(layerHeight)
        hull(){
            translate([radius, 0])
            circle(cornerRadius);
            translate([radius - cornerLength, 0])
            circle(cornerRadius);
        }
}


module layer(radius, height) {
        outerPoints = getPoints(angles, radius);
        innerPoints = getPoints(angles, radius-2);
        linear_extrude(height, twist = twist)
        difference(){
            polygon(points = outerPoints);
            polygon(points = innerPoints);
        }
}

function getPoints(angles, radius) = [for (th=angles) [radius*cos(th), radius*sin(th)]];
function getAngleAdjustment(steps) = -(steps * angleIncrement) % 360;
function getSinRadius(currentSegment, interval) = fixedRadius + sinRadius * sin(180*(currentSegment/segmentPatternSteps[interval]));
function getNestedRadius(currentSegment, interval) = fixedRadius * (sizeRatio ^ (currentSegment * 2.5));
function angleFromTriangleSides(op, aj1, aj2) = acos ((aj1 ^ 2 + aj2 ^ 2 - op ^ 2) / ( 2 * aj1 * aj2));

module stackedPolygons(){
    for (corner = [0:sides]){
        rotate(corner * polygonAngle)
        for (i=[0:hullInterval:(segmentSteps[0] - hullInterval)]){
            hull(){
                for (j=[i:i+hullInterval]){
                    variableRatio = - cos(180*(j/segmentPatternSteps[0]));
                    radiusOffset = - (cornerFixedInset + cornerVariableInset * variableRatio);
                    angleOffset = cornerFixedOffsetAngle + variableRatio * cornerVariableOffsetAngle;
                    rotate([0, 0, getAngleAdjustment(j) + angleOffset])
                    translate([0, 0, j*layerHeight])
                        corner(getSinRadius(j, 0) + radiusOffset, layerHeight);
                }
            }
        }
        rotate(corner * polygonAngle)
        for (i=[segmentEndSteps[0]:hullInterval:segmentEndSteps[1]]){
            hull(){
                for (j=[i:i+hullInterval]){
                    variableRatio = (1 - (j -segmentEndSteps[0])/segmentPatternSteps[1]);
                    radiusOffset = - (cornerFixedInset + cornerVariableInset * variableRatio);
                    angleOffset = cornerFixedOffsetAngle + variableRatio * cornerVariableOffsetAngle;
                    rotate([0, 0, getAngleAdjustment(j)+ angleOffset])
                    translate([0, 0, j*layerHeight])
                        corner(getNestedRadius(j - segmentEndSteps[0], 1) + radiusOffset ,layerHeight);
                }
            }
        }
    }
    layerInterval = 2;
    for (j=[0:segmentSteps[0]]){
        if (j % layerInterval == 0){
            rotate([0, 0, getAngleAdjustment(j)])
            translate([0, 0, j*layerHeight])
                layer(getSinRadius(j, 0), layerHeight);
        }
    }
    for (j=[segmentEndSteps[0] : segmentEndSteps[1]-70]){
        if (j % layerInterval == 0){
            rotate([0, 0, getAngleAdjustment(j)])
            translate([0, 0, j*layerHeight])
                layer(getNestedRadius(j-segmentSteps[0], 1), layerHeight);
        }
    }
}


stackedPolygons();