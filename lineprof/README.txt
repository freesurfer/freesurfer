
liblineprof 1.0

Martin Reuter, Oct 17 2012

This library solves the Laplace equation in a 2D piecewise linear
polygon with four boundary segments:

0-level, 1-level and the two connecting boundaries (segmentL, segmentR):

 1111111111111111111
 L                 R
 L                 R
 0000000000000000000

It builds against ITK, VTK and Petsc.

The library extends code by Gheorghe Postelnicu, 2006 (pyScout).



Description:
============

1. The user creates and passes the polygon with 4 boundary segments as 
an array of 2D points:

   std::vector < std::vector < double > > points2d

where the inner vector's size is 2 for x and y coordianate, 
and 4 integer arrays for the boundary segments:

  std::vector < int > segment0
  (similar segment1, segmentLeft, segmentRight)
  
with indices into the points2d array (starting with zero). Both the 
segment0 and segment1 are indexed from left to right, the segmentLeft
and segmentRight from the 0-level to the 1-level.


2. The library then solves the Laplace equation with appropriate boundary
contitions using an image as grid to construct the Laplace matrix and right
handside vector and then solves the problem using Petsc.


3. It then traces line profiles from the 0-level to the 1-level at the 
requested density and returns those lines as an array, with each
line as an array of 2D point coordinates:

  std::vector < std::vector < std::vector < double > > > profiles

where the inner vector's size is 2 for x and y coordinate.



See LineProfTest.cpp for an example.

