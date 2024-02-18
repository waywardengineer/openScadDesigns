
sides = 3;
layerHeight = 0.2;
twistPer100mm = 100;
baseRadius = 25; //min 24
outerRadius = 35; // min 26

polygonAngle = 360 / sides;
angleIncrement = twistPer100mm * (layerHeight / 100);


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
function zeroNegatives(list) = [for (item = list) max(0, item)];
    
module circularScrewSlot(headDia, slotDia, angleSpan, radius){
    for (angle=[0:0.1:angleSpan]){
        rotate(angle)
        translate([radius, 0, 0])
        circle(slotDia / 2);
    }
    translate([radius, 0, 0])
    circle(headDia / 2);
}

module bottomCap(radius){
    holeRadius = 20;
    coneHeight = 2;
    wallThickness = 0.1;
    screwThreadDia = 2.7;
    screwHeadDia = 5;
    screwSlotAngleSpan = 20;
    screwSlotCircleDia = 30;
    difference(){
        union(){
            linear_extrude(1.2)
            difference(){
                layerP_CR(radius, radius, 0, 0);
                rotate(-7)
                union(){
                    for (slotAngleInterval = [0:3]){
                        rotate(slotAngleInterval * 120)
                        circularScrewSlot(screwHeadDia, screwThreadDia, screwSlotAngleSpan, screwSlotCircleDia/2);
                    }
                }
            }
            hull(){
                translate([0, 0, coneHeight])
                cylinder(0.01, holeRadius/2 + wallThickness * sqrt(2));
                cylinder(0.01, holeRadius/2 + wallThickness * sqrt(2) + coneHeight);            
            }        
        }
        hull(){
            translate([0, 0, coneHeight+0.01])
            cylinder(0.01, holeRadius/2);
            translate([0, 0, -0.01])
            cylinder(0.01, holeRadius/2 + coneHeight);            
        }
    }

}
      



module spherical(){
    sphereRadius = outerRadius;
    heightOffset = sphereRadius - sqrt(sphereRadius ^ 2 - baseRadius ^ 2);
    R1 = [for (height = [heightOffset:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - (sphereRadius-height) ^2)];
    R2 = [for (height = [0:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - height ^ 2)];
    R = cat(R1, R2);
    nR = len(R);
    int = [0:nR-1];
    
    /*PRt1 = nR-21;
    PRint1 = [0:PRt1-1];
    PRint2 = [PRt1: nR-1];
    PR1 = [for (i = PRint1) (i % 7 < 3) ? R[i] : 0]; 
    PR2 = [for (i = PRint2) (i % 5 < 3) ? R[i] : 0]; 
    PR = cat(PR1, PR2);*/
    PR = [for (i = int) (i % 1 < 2) ? R[i] : 0]; 
        
    
    PWt1 = nR-21;
    PWint1 = [0:PWt1-1];
    PWint2 = [PWt1: nR-1];
    PW1 = [for (i = PWint1) 1]; 
    PW2 = [for (i = PWint2) 1 + (i-PWt1) / 10]; 
    PW = cat(PW1, PW2);
    //PW = [for (i = int) 1];
        
    RO = [for (i = int) R[i]-0.5];
    
    RIt1 = round(len(R) * 0.5);
    RIint1 = [0:RIt1-1];
    RIint2 = [RIt1: nR-1];
    RI1 = [for (i = RIint1) RO[i] - 2];
    RI2 = [for (i = RIint2) RO[i] - (2 + (i-RIt1) * layerHeight * 0.7)];
    RI = cat(RI1, zeroNegatives(RI2));
        
    LR = [for (i = int) (i * angleIncrement) % 360];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        layerP_CR(PR[i], PW[i], RO[i], RI[i], 0.8);
    }
}

module sinLongTip(){
    
    R1 = [for (height = [0:layerHeight:60])  baseRadius + (outerRadius - baseRadius) * sin(180*(height/60))];
    R2 = [for (height = [0:layerHeight:100])  baseRadius * (0.4 ^ ((height / 100) * 2.5))];
    R = cat(R1, R2);
    nR = len(R);
    int = [0:nR-1];
    
    /*PRt1 = nR-21;
    PRint1 = [0:PRt1-1];
    PRint2 = [PRt1: nR-1];
    PR1 = [for (i = PRint1) (i % 7 < 3) ? R[i] : 0]; 
    PR2 = [for (i = PRint2) (i % 5 < 3) ? R[i] : 0]; 
    PR = cat(PR1, PR2); 
    
    PW = [for (i = int) 1];*/
        
    RO = [for (i = int) R[i]];
    
    RI = zeroNegatives([for (i = int) RO[i] - 2]);
        
    LR = [for (i = int) (360/(2.3*sides)) * sin(180*(i/nR))];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        layer_R(sides, RO[i], RI[i], 1, 1);
    }
    for (i = int){
        rotate([0, 0, -LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        layer_R(sides, RO[i], RI[i], 1, 1);
    }
}


bottomCap(baseRadius+1);
spherical();