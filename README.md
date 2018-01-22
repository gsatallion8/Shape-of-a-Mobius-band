# Shape-of-a-Mobius-band
A solver to compute the shape of a mobius band for different material properties and aspect ratios. 
The Kirchoff-Love equations (12th order ODE) governing the shape of slender elastic members are solved numerically using Euler explicit scheme. 


In the non-dimensional regime, the two parameters \alpha and \beta (a and b in the code) are the only two parameters that define the problem. Hence, the code generates the shape (centerline coordinates and inclination of the crosssection) for various values of these two parameters.
