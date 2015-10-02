#include "declarations.h"

double CircleOverlap(double d, double r1, double r2)
{
  double r, R, Area, r_2, R_2, d_2;
  
  //d = sqrt(sqr(x2-x1) + sqr(y2-y1));
  
  if(r1 > r2){
    R = r1;
    r = r2;
  }
  else{
    R = r2;
    r = r1;    
  }
  
  if(d >= (R+r)) Area = 0.0;
  else if(d <= (R-r))  Area = PIE * sqr(r);
  else{
    r_2 = sqr(r);
    R_2 = sqr(R);
    d_2 = sqr(d);
    
    Area = (r_2*acos( (d_2 + r_2 - R_2) / (2.0 * d * r) ) +
	    R_2*acos( (d_2 + R_2 - r_2) / (2.0 * d * R) ) -
	    0.5 * sqrt((r + R - d) * (r - R + d) * (R - r + d) * (r + R + d)));
  }
  
  return(Area);
}
