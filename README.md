# QuadraticProgramming (Under development)
Created by Dominic Keehan [Majority of initial development in May 2020]

This repository contains a complete implemenation of a linear program solver, including both simplex and revised simplex method solvers ("fullsimplex.m" and "RSM.m" respectively). Running documentation can be found in the functions themselves. These solvers were developed for a university course in operations research and allow a user to solve the problem "minimise cx s.t. Ax=b, x>=0". The goal of this was to develop a solid understanding of the simplex and revised simplex method and linear programming as a whole (convexity, etc...). However, this course did not touch on quadratic programming, solving problems of the form "minimise cx + 1/2x'Cx, s.t. Ax=b, x>=0". A new goal is to gain an understanding of how these problems are solved, using "Wolfe's method" and extend the current linear solvers to meet this functionality.

http://pages.cs.wisc.edu/~brecht/cs838docs/wolfe-qp.pdf

Goals:

Develop an understanding of the intricacies of quadratic programming

Future work:

Get the Quadratic solver running
