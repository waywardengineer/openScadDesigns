
sides = 12;
layerHeight = 0.2;
twistPer100mm = 100;
baseRadius = 30;
outerRadius = 40;

polygonAngle = 360 / sides;
angleIncrement = twistPer100mm * (layerHeight / 100);


module flatLayer(PR, PW, RO, RI, RR) {
    if (PR > 0){
        difference(){
            flatPolygon(PR);
            flatPolygon(PR - PW);          
        }
    }
    if (RO > 0){
        for (corner = [0:sides]){
            rotate(corner * polygonAngle)
            hull(){
                translate([RO, 0])
                circle(RR);
                translate([RI, 0])
                circle(RR);
            }
        }
        
    }
}

module flatPolygon(radius){
    union(){
        polygon(points = getPoints(getAngles(sides), radius));
        //polygon(points = getPoints(getAngles(sides*3), radius*0.95));
    }
}

function getAngles(sides) = [for (j = [0:sides-1]) j*(360/sides)];
function getPoints(angles, radius) = [for (th=angles) [radius*cos(th), radius*sin(th)]];
function angleFromTriangleSides(op, aj1, aj2) = acos ((aj1 ^ 2 + aj2 ^ 2 - op ^ 2) / ( 2 * aj1 * aj2));
function cat(L1, L2) = [for(L=[L1, L2], a=L) a];
function zeroNegatives(list) = [for (item = list) max(0, item)];

module bottomCap(radius){
    holeRadius = 25;
    coneHeight = 2;
    wallThickness = 1.2;
    screwThreadDia = 2.7;
    screwHeadDia = 5;
    screwSlotAngleSpan = 20;
    screwSlotCircleDia = 40;
     difference(){
        union(){
            linear_extrude(wallThickness)
                difference(){
                    flatLayer(radius, radius, 0, 0);
                    for (slotAngleInterval = [0:3]){
                        rotate(slotAngleInterval * 120)
                        union(){
                            for (angle=[0:0.1:screwSlotAngleSpan]){
                                rotate(angle)
                                translate([screwSlotCircleDia/2, 0, 0])
                                circle(screwThreadDia / 2);
                            }
                            translate([screwSlotCircleDia/2, 0, 0])
                            circle(screwHeadDia / 2);
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
    
    PRt1 = nR-21;
    PRint1 = [0:PRt1-1];
    PRint2 = [PRt1: nR-1];
    PR1 = [for (i = PRint1) (i % 7 < 3) ? R[i] : 0]; 
    PR2 = [for (i = PRint2) (i % 5 < 3) ? R[i] : 0]; 
    PR = cat(PR1, PR2); 
    
    PW = [for (i = int) 1];
        
    RO = [for (i = int) R[i] - 0.75];
    
    RIt1 = round(len(R) * 0.7);
    RIint1 = [0:RIt1-1];
    RIint2 = [RIt1: nR-1];
    RI1 = [for (i = RIint1) RO[i] - 2];
    RI2 = [for (i = RIint2) RO[i] - (2 + (i-RIt1) * layerHeight * 0.8)];
    RI = cat(RI1, zeroNegatives(RI2));
        
    LR = [for (i = int) (i * angleIncrement) % 360];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        flatLayer(PR[i], PW[i], RO[i], RI[i], 0.8);
    }
}

module sinLongTip(){
    
    R1 = [for (height = [0:layerHeight:60])  baseRadius + (outerRadius - baseRadius) * sin(180*(height/60))];
    R2 = [for (height = [0:layerHeight:100])  baseRadius * (0.4 ^ ((height / 100) * 2.5))];
    R = cat(R1, R2);
    nR = len(R);
    int = [0:nR-1];
    
    PRt1 = nR-21;
    PRint1 = [0:PRt1-1];
    PRint2 = [PRt1: nR-1];
    PR1 = [for (i = PRint1) (i % 7 < 3) ? R[i] : 0]; 
    PR2 = [for (i = PRint2) (i % 5 < 3) ? R[i] : 0]; 
    PR = cat(PR1, PR2); 
    
    PW = [for (i = int) 1];
        
    RO = [for (i = int) R[i] - 0.75];
    
    RI = zeroNegatives([for (i = int) RO[i] - 2]);
        
    LR = [for (i = int) (i * angleIncrement) % 360];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        flatLayer(PR[i], PW[i], RO[i], RI[i], 0.8);
    }
}


bottomCap(baseRadius);
sinLongTip();