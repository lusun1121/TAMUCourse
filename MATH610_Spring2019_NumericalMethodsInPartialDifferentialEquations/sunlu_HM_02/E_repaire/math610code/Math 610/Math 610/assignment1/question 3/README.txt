
This is a 1D Finite Element Code written for the Texas A&M Math 610 Lab
course on numerical analysis of partial differential equations (FEM)
in a manner that shows how a finite element code could look when written 
with the goal of modularity, flexibility and design structure.  It is not
perfect nor highly optimized but it does separate the various tasks into
their own modules and attempts to minimize the reuse of code.  There are
undoubtedly many ways to speed up and reduce the memory footprint, but this
is more a work for teaching structure and helping the student see how to
break the tasks down into simple pieces that fit together to solve the full
system.  

In 1D, the code is typically so simple that building in the proper 
structure seems unnecessary. However, the jump from 1D to 2D is often 
problematic and somewhat more challenging than it should be because we
don't learn good habits in the first place.

To this end, we have written this code with the habits that a higher 
dimensional code would employ to be successful. This includes proper
comments, breaking the problem into smaller sub pieces which allow you to
focus on a certain part of the entire piece at a time.  Also, we have 
used some more advanced MATLAB tools called structured arrays (type 'help
struct' in command window to see more about this) to group different parts 
of the code together.  We do this for the triangulation, the Quadrature on
the reference elements, the dof handler, and the finite element bases
in the reference element evaluated at the corresponding quadrature points
on the reference element.

This code implements the laplacian with Dirichlet boundary 
conditions on the interval Lambda = (0,L). 

The boundary condition locations are set in the constructTriangulation.m
file.  The Quadrature on the bulk reference elements are set in 
getQuadOnRefElement.m  You may have to write your own Quadrature on edge
reference elements code if your problem requires integrals over edges of 
elements. The degrees of freedom are constructed from the triangulation 
and stored in the constructDoFHandler.m function.  The meat of the work is
 done in the assembleLocal*.m files which  construct the matrix and rhs 
terms on each local element based on the finite elements chosen.  

The whole system is put together and run from the mainFEM1D_Dirichlet.m 
script file which can be executed.  We have also turned that mainscript
into a function called runFEM1d.m which takes in n_elements, p polynomial
basis degree, quad_degree and a flag to suppress output of plots and 
returns the errors in the various norms.  Then we provide a script called 
mainErrorAnalysis.m which does an error analysis of code by calling the 
function with multiple times with increasing number of elements (doubling
each time) and thus, decreasing mesh size.  The one step rate of 
convergence is computed and output in a table format with the errors in 
various norms and the corresponding rates.

Good luck playing around with this code.  I hope it is useful to learn good 
structure and basic design of fem codes.

Spencer Patty
Jan 24, 2017



