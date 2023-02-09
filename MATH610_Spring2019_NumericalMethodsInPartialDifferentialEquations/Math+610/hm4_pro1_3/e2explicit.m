function e2explicit(tau_size)%foward
L=3;
%tau_size=19;%22;
%tau_size=40;%43;
%tau_size=79;%82;
tau=L/(tau_size-1);
t=linspace(0,L,tau_size);
fprintf('tau=%d,tau size=%d\n',tau,tau_size);


Q=1;
p = 1;  % polynomial degree of lagrange finite element basis
%u_exact = @(x,y,t)t.*exp(-t).*cos(3*pi*x).*cos(pi*y);%nx1
%grad_u_exact = @(x,y,t) [-3*pi*t.*exp(-t).*sin(3*pi*x).*cos(pi*y),...
%                         -pi*t.*exp(-t).*cos(3*pi*x).*sin(pi*y)];%nx2
% right hand side function
f = @(x,y,t)ones(size(x));
%f = @(x,y,t)(1-t+10*pi^2*t).*exp(-t).*cos(3*pi*x).*cos(pi*y)+...
%    Q*u_exact(x,y,t);%nx1 
%g_D = @(x,y,t)0;%@(x,y,t) u_exact(x,y,t);

quad_n_points = 5;
Quad = getQuadOnRefElement(quad_n_points);
[FE_at_Quad] = feEval( Quad, p );

%get triangle and degree of freedom
%output=zeros(3,4);
for j=1:3
    switch j
        case 1
            m_size = 10;
        case 2
            m_size = 20;
        case 3
            m_size = 40;
    end
    
    T=gettriangle(j);

    DoFHandler = constructDoFHandler(T,p);

    M = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes);
    S = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes);

    % assemble_system
    for cell = 1:T.n_elements
        dofIndices = DoFHandler.dofs(cell,:); % [1x(p+1)]  extract indices pertaining to cell nodes
        vertices = T.nodes(T.elements(cell,:),:); % [(dim+1)xdim] coordinates of vertices
        [mass_matrix,stiffness_matrix] = local_assembly(vertices, FE_at_Quad, Quad, p); 
    
        % contribute local terms to global stiffness and RHS structures
        M(dofIndices,dofIndices) = M(dofIndices,dofIndices) + mass_matrix;
        S(dofIndices,dofIndices) = S(dofIndices,dofIndices) + stiffness_matrix;
    end

    uh1 = zeros(DoFHandler.n_dofs,1);
    %uh2 =  zeros(DoFHandler.n_dofs,1);


    %forward Euler
    for i=1:tau_size
        ti=t(i);
        RHS = zeros(DoFHandler.n_dofs,1);
        for cell = 1:T.n_elements
            dofIndices = DoFHandler.dofs(cell,:); 
            vertices = T.nodes(T.elements(cell,:),:); 
            local_rhs=local_assemblyrhs(vertices,f,ti,FE_at_Quad, Quad, p);
            RHS(dofIndices)=RHS(dofIndices)+local_rhs;
        end
        uh2=M\(((1-tau*Q)*M-tau*S)*uh1+tau*RHS);
        %uh2 = ((1-tau+Q*tau)*M+tau*S)\(M*uh1+tau*RHS);
        %uh2=(1-tau)*uh1-tau*(M\S)*uh1+tau*M\RHS;
        %uh2=(1-tau)*uh1-tau*(M\S)*uh1+tau*M\RHS;

        %uh2=(1-Q*tau)*uh1-tau*(M\S)*uh1+tau*M\RHS;
        %compute error
        if (ti==1 || ti==3)
            lmbd=figure;
            plotfunction(uh2,T,DoFHandler); 
            title(['m=',num2str(m_size),',t=',num2str(ti),',tau=',num2str(tau),',tau size=',num2str(tau_size)]);
            %saveas(lmbd, ['C:\Users\Sunlu\Desktop\HM_610\hm4_pro1_3\hm4_prob2_ex_m',...
            %    num2str(m_size),'_t',num2str(ti),'_tau',num2str(tau_size),'.png'],'png');

            %if (i==5 || i==10)
            %Quad_Error = getQuadOnRefElement(quad_n_points);
            %FE_at_Quad_Error = feEval(Quad_Error, p);
        
            %L2sqrd = 0;
            %H1sqrd = 0;
            %for cell = 1:T.n_elements
            %    dofIndices = DoFHandler.dofs(cell,:); 
            %    vertices = T.nodes(T.elements(cell,:),:); 
    
             %   [localL2sqrd, localH1sqrd] = computeLocalErrors(vertices, uh1(dofIndices),...
             %       u_exact, grad_u_exact,ti,Quad_Error, FE_at_Quad_Error, p);
             %   L2sqrd = L2sqrd + localL2sqrd;
             %   H1sqrd = H1sqrd + localH1sqrd;
            %end
            %L2 = sqrt(L2sqrd);
            %H1 = sqrt(H1sqrd);
            %if ti==1
            %    output(j,1)=L2;%t=1, L2 error
            %    output(j,2)=H1;%t=1, H1 error
            %else
            %    output(j,3)=L2;%t=3, L2 error
            %    output(j,4)=H1;%t=3, H1 error
            %end
        end
        uh1=uh2;
    end
end
%h=[1/10^2,1/20^2,1/40^2];

%for pt = 1:3
%    fprintf('%f & %f & %f & %f & %f\\\\\n',h(pt),output(pt,1),output(pt,2),output(pt,3),output(pt,4));
%end
end