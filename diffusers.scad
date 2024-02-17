
sides = 12;
layerHeight = 0.2;
//fixedRadius = 30;
//sinRadius = 60;
twistPer100mm = 100;
twist = 0;
bottomRadius = 30;
sphereRadius = 40;

cornerRadius = 0.8;

polygonAngle = 360 / sides;
angleIncrement = twistPer100mm * (layerHeight / 100);
sideIntersectionRatio =  sin(polygonAngle - angleIncrement) / sin(angleIncrement);
sizeRatio = sin(180 - polygonAngle) / ((1 + sideIntersectionRatio) * sin(angleIncrement));


module flatLayer(PR, PW, RO, RI) {
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
                circle(cornerRadius);
                translate([RI, 0])
                circle(cornerRadius);
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


module stackedPolygons(){
    
    heightOffset = sphereRadius - sqrt(sphereRadius ^ 2 - bottomRadius ^ 2);
    R1 = [for (height = [heightOffset:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - (sphereRadius-height) ^2)];
    R2 = [for (height = [0:layerHeight:sphereRadius])  sqrt(sphereRadius ^ 2 - height ^ 2)];
    R = cat(R1, R2);
    nR = len(R);
    int = [0:nR-1];
    
    t1 = nR-21;
    int1 = [0:t1-1];
    int2 = [t1: nR-1];
    PR1 = [for (i = int1) (i % 7 < 3) ? R[i] : 0]; 
    PR2 = [for (i = int2) (i % 5 < 3) ? R[i] : 0]; 
    PR = cat(PR1, PR2); 
    
    PW = [for (i = int) 1];
        
    RO = [for (i = int) R[i] - 0.75];
    
    t2 = round(len(R) * 0.7);
    int3 = [0:t2-1];
    int4 = [t2: nR-1];
    RI1 = [for (i = int3) RO[i] - 2];
    RI2 = [for (i = int4) RO[i] - (2 + (i-t2) * layerHeight * 0.8)];
    RI = cat(RI1, zeroNegatives(RI2));
        
    LR = [for (i = int) (i * angleIncrement) % 360];

    for (i = int){
        rotate([0, 0, LR[i]])
        translate([0, 0, i*layerHeight])
        linear_extrude(layerHeight)
        flatLayer(PR[i], PW[i], RO[i], RI[i]);
    }
}


bottomCap(bottomRadius);
stackedPolygons();