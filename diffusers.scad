
sides = 3;
layerHeight = 0.2;
twistPer100mm = 100;

polygonAngle = 360 / sides;
angleIncrement = twistPer100mm * (layerHeight / 100);
baseRadii = [24, 25, 25, 25, 25];
sphereRadii = [30, 35, 40, 45, 50];
hSpacing = 0.3;
vSpacing = 0.2;
cornerRadius = 0.6;
bulbRadius = 10;
bottomThickness = 1.2;


module layerP_CR(PR, PW, RO, RI, CR) {
    if (PR > 0){
        difference(){
            union(){
                layer_P(PR, sides);
                layer_P(PR * 0.55, sides * 2);            
            }
            union(){
                layer_P(max(PR - PW, 0), sides);
                layer_P(max(0, PR * 0.55 - PW), sides * 2);
            }
        }
    }
    if (RO > 0){
        layer_R(sides, RO, RI, CR, CR);
    }
}


module layer_R(sides, RO, RI, CRO, CRI){
    for (angle = getAngles(sides)){
        rotate(angle)
        hull(){
            translate([RO, 0])
            circle(CRO);
            translate([RI, 0])
            circle(CRI);
        }
    }
}

module layer_P(radius, sides){
    union(){
        polygon(points = getPoints(getAngles(sides), radius));
    }
}

function getAngles(sides) = [for (j = [0:sides-1]) j*(360/sides)];
function getPoints(angles, radius) = [for (th=angles) [radius*cos(th), radius*sin(th)]];
function angleFromTriangleSides(op, aj1, aj2) = acos ((aj1 ^ 2 + aj2 ^ 2 - op ^ 2) / ( 2 * aj1 * aj2));
function cat(L1, L2) = [for(L=[L1, L2], a=L) a];
function cat3(L1, L2, L3) = cat(cat(L1, L2), L3);
function zeroNegatives(list) = [for (item = list) max(0, item)];
    

module circularScrewHoleFlatPattern(patternDia, headDia, threadDia, numHoles, slotAngleSpan, patternRotation){
    for (slotAngleInterval = [0:numHoles-1]){
        rotate(slotAngleInterval * 120 + patternRotation)
        union(){
            for (angle=[0:0.1:slotAngleSpan]){
                rotate(angle)
                translate([patternDia / 2, 0, 0])
                circle(threadDia / 2);
            }
            translate([patternDia / 2, 0, 0])
            circle(headDia / 2);
        }
    }
}
      
module bulbOpening(bodyHeight){
    openingHeight = bodyHeight - 15;
    h1 = 10;
    h1_2 = bulbRadius * 0.5;
    h2_3 = bulbRadius * 1.5;
    h2 = bodyHeight - h1 - h1_2 - h2_3;
    
    
    r1 = bulbRadius;
    r2 = bulbRadius * 1.5;
    
    difference(){
        union(){
            cylinder(h1, r1, r1);
            translate([0, 0, h1])
            cylinder(h1_2, r1, r2);
            translate([0, 0, h1 + h1_2])
            cylinder(h2, r2, r1);
            translate([0, 0, h1 + h1_2 + h2])
            cylinder(h2_3, r1, 0);

        }
        cylinder(bodyHeight, 2, 2);
    }    
}

module bulb(sizeIndex, parts){
    sphereRadius = sphereRadii[sizeIndex];
    baseRadius = baseRadii[sizeIndex];
    heightOffset = sphereRadius - sqrt(sphereRadius ^ 2 - baseRadius ^ 2);
    R2 = [for (height = [heightOffset:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - (sphereRadius-height) ^2)];
    R3 = [for (height = [0:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - height ^ 2)];
    R = cat(R2, R3);
    nR = len(R);
    int = [0:nR-1];
    if (parts < 2){
        difference(){
            ribs(R, nR, int, 0, 0);
            cylinder(bottomThickness, bulbRadius, bulbRadius);
            linear_extrude(bottomThickness)
            circularScrewHoleFlatPattern(30, 5, 2.7, 3, 20, -7);
            translate([0, 0, bottomThickness])
            linear_extrude(8)
            circularScrewHoleFlatPattern(30, 5, 5, 3, 20, -7);

        }
    }
    if (parts == 0 || parts == 2){
        difference(){
            body(R, nR, int);
            ribs(R, nR, int, 0.1, 0.1); 
            bulbOpening(nR * layerHeight);
            translate([0, 0, bottomThickness])
            linear_extrude(8)
            circularScrewHoleFlatPattern(30, 5, 5, 3, 20, -7);
        }
    }
}

module ribs(R, nR, int, hSpacing, vSpacing){    
    PRt1 = bottomThickness/layerHeight;
    PRint1 = [0:PRt1-1];
    PRint2 = [PRt1: nR-1];
    PR1 = [for (i = PRint1) R[i] + 3 * hSpacing]; 
    PR2 = [for (i = PRint2) 0]; 
    PR = cat(PR1, PR2);
    
    PW = [for (i = int) PR[i] - bulbRadius];
    
    RO = [for (i = int)  R[i] -0.5];
        
    RI = zeroNegatives([for (i = int) RO[i] - (10 + 3 * hSpacing - 2 * (i/nR))]);
        
    LR = [for (i = int) (i * angleIncrement) % 360];
    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight - vSpacing])
        linear_extrude(layerHeight + vSpacing * 2)
        layerP_CR(PR[i], PW[i], RO[i], RI[i],  cornerRadius + hSpacing);
    }
}


module body(R, nR, int){        
    
    PW = [for (i = int) R[i]];
        
    RO = [for (i = int) 0];
    
    RI = [for (i = int) 0];
        
    LR = [for (i = int) (i * angleIncrement) % 360];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        layerP_CR(R[i], PW[i], RO[i], RI[i], 0.4);
    }
}
    //bulb(4);

    bulb(4, 2);

