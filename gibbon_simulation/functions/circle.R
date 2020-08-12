# functions to compute properties of overlapping circles. 
# (x, y, r) defines the circle at point (x, y) of radius r.

# function to compute intersection of circle_1 and circle_2 and identify which lies in circle 3
circle_intersect = function(x1, y1, r1, x2, y2, r2, x3, y3, r3){
  d = ((x2 - x1)^2 + (y2 - y1)^2)^0.5
  a = (r1^2 - r2^2 + d^2) / (2 * d)
  h = (abs(r1^2 - a^2))^0.5
  xx = x1 + (a * (x2 - x1)) / d
  yy = y1 + (a * (y2 - y1)) / d
  x12 = xx + (h * (y2 - y1) / d)
  y12 = yy - (h * (x2 - x1) / d)
  if((x12 - x3)^2 + (y12 - y3)^2 > r3^2){
    x21 = xx - (h * (y2 - y1) / d)
    y21 = yy + (h * (x2 - x1) / d)
    return(c(x21, y21))
  } else{
    return(c(x12, y12))
  }
}

# function to compute magnitude of area common to three circles (x1, y1, r1),...,(x3, y3, r3).
common_area = function(x1, y1, r1, x2, y2, r2, x3, y3, r3){
  int12 = circle_intersect(x1, y1, r1, x2, y2, r2, x3, y3, r3)
  int23 = circle_intersect(x2, y2, r2, x3, y3, r3, x1, y1, r1)
  int31 = circle_intersect(x3, y3, r3, x1, y1, r1, x2, y2, r2)
  a = ((int31[1] - int12[1])^2 + (int31[2] - int12[2])^2)^0.5
  b = ((int23[1] - int12[1])^2 + (int23[2] - int12[2])^2)^0.5
  c = ((int31[1] - int23[1])^2 + (int31[2] - int23[2])^2)^0.5
  s = 0.5 * (a + b + c)
  A = c(a, b, c)
  R = c(r1, r2, r3)
  return(sum(R^2 * asin(A / (2*R)) - (A / 4)*(4*R^2 - A^2)^0.5) + (s * (s - a) * (s - b) * (s - c))^0.5)
}

# function to compute area common to the two circles (x1, y1, r1), (x2, y2, r2)
intersect_area = function(x1, y1, r1, x2, y2, r2){
  d = ((x2 - x1)^2 + (y2 - y1)^2)^0.5
  return((r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1)) + r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))) -
           0.5 * ((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))^0.5)
}
