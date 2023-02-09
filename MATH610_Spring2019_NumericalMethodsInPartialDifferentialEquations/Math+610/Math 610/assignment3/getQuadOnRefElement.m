function [Quad] = getQuadOnRefElement(num_q_points)
%xhat  nq*2
%what  nq*1
switch(num_q_points)

    case 1,
        Quad.nq = 1;

        Quad.xhat = [0.333333333333333,0.333333333333333];

        Quad.what = [0.500000000000000];

    case 2,
    
        Quad.nq = 3;

        Quad.xhat =[0.166666666666667,0.166666666666667;...
                    0.166666666666667,0.666666666666667;...
                    0.666666666666667,0.166666666666667]; 
    
        Quad.what =[0.166666666666667;...
                    0.166666666666667;...
                    0.166666666666667];

    case 3,
      
        Quad.nq = 4;
   
        Quad.xhat = [ 0.333333333333333,0.333333333333333;... 
                      0.200000000000000,0.200000000000000;...
                      0.200000000000000,0.600000000000000;...
                      0.600000000000000,0.200000000000000];
      
        Quad.what = [ -0.281250000000000;...
                      0.260416666666667;...
                      0.260416666666667;...
                      0.260416666666667];

    case 4,

        Quad.nq = 6;
   
        Quad.xhat = [0.091576213509771,0.091576213509771;...
                     0.091576213509771,0.816847572980459;...
                     0.816847572980459,0.091576213509771;...
                     0.445948490915965,0.445948490915965;...
                     0.108103018168070,0.445948490915965;...
                     0.445948490915965,0.445948490915965];
 
        Quad.what = [0.054975871827661;...
                     0.054975871827661;...
                     0.054975871827661;...
                     0.111690794839005;...
                     0.111690794839005;...
                     0.111690794839005;];
    case 5,

        Quad.nq = 7;

        Quad.xhat = [0.333333333333333,0.333333333333333;... 
                     0.797426985353087,0.101286507323456;...
                     0.101286507323456,0.797426985353087;...
                     0.101286507323456,0.101286507323456;...
                     0.059715871789770,0.470142064105115;...
                     0.470142064105115,0.059715871789770;...
                     0.470142064105115,0.470142064105115];

        Quad.what = [0.112500000000000;...
                     0.062969590272414;...
                     0.062969590272414;...
                     0.062969590272414;...
                     0.066197076394253;...
                     0.066197076394253;...
                     0.066197076394253];
   
    otherwise,
        Quad.nq = 0;
        Quad.what = [];
        Quad.xhat = [];
end




