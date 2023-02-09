
This is a 1D Finite Element Code written for the Texas A&M Math 610 Lab
course on numerical analysis of partial differential equations (FEM)
in a manner that shows how a dimensional code would look when written well.
In 1D, the code is typically so simple that building in the proper 
structure seems unnecessary. However, the jump from 1D to 2D is often 
problematic and somewhat more challenging that it should be because we
don't learn proper habits in the first place.

To this end, we have written this code with the habits that a higher 
dimensional code would employ to be successful. This includes proper
comments, breaking the problem into smaller sub pieces which allow you to
focus on a certain part of the entire piece at a time.  Also, we have 
used some more advanced MATLAB tools called structured arrays (type 'help
struct' in command window to see more about this) to group different parts 
of the code together.  We do this for the triangulation, the Quadrature on
the reference elements, and the finite elements evaluated at the quadrature
points on the reference elements.

This code implements the laplacian with Dirichlet boundary 
conditions on the interval Lambda = (0,L). 

The boundary condition locations are set in the constructTriangulation.m
file.  The Quadrature on the bulk reference elements are set in 
getQuadOnRefElement.m  YOu may have to write your own Quadrature on edge
reference elements code if your problem requires integrals over edges of 
elements.  The meat of the work is done in the assembleLocal*.m files which 
construct the matrix and rhs terms on each local element based on the
finite elements chosen.  

The whole system is put together and run from the mainFEM1D_Dirichlet.m 
script file which can be executed or can be rewritten into a function that
takes in n, the number of elements, and returns the errors for that given 
n.  This might be easier when you are running an error analysis to verify 
the convergence rates of the FEM in the different norms.  If you go that 
route, you may also want to add a flag to output plots of solution and 
errors or not since when running error analysis, the plots sometimes get 
in the way. It is up to you!  

Good luck playing around with this code.  I hope it is useful to learn good 
structure and basic design of fem codes.

Spencer Patty
Jan 24, 2017



