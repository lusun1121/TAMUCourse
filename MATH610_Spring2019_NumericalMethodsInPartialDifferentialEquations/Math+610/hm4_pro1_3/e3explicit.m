function e3explicit(tau_size)%foward
L=3;
%tau_size=19;%22;
%tau_size=40;%43;
%tau_size=79;%82;
tau=L/(tau_size-1);
t=linspace(0,L,tau_size);
fprintf('tau=%d,tau size=%d\n',tau,tau_size);


Q=1;
p = 1;  
% right hand side function
f = @(x,y,t)x.*y;%nx1 
g= @(x)1;


quad_n_points = 5;
Quad = getQuadOnRefElement(quad_n_points);
[FE_at_Quad] = feEval( Quad, p );

Quad1=getQuadOnRefElement1D(quad_n_points);
[FE_at_Quad1] = feEval1D( Quad1, p );


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
    %uh2 = zeros(DoFHandler.n_dofs,1);


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
        
        for cell=1:T.n_boundaryedges
            dofIndices=T.edges(T.boundaryedges(cell),:);
            v1=T.nodes(T.edges(T.boundaryedges(cell),1),:);
            v2=T.nodes(T.edges(T.boundaryedges(cell),2),:);
            [rhs_g] = assembleLocalRhsg(g,v1,v2, FE_at_Quad1, Quad1);
            RHS(dofIndices) = RHS(dofIndices) +rhs_g;
        end
        
        uh2=M\(((1-tau*Q)*M-tau*S)*uh1+tau*RHS);
        %uh2 = ((1-tau+Q*tau)*M+tau*S)\(M*uh1+tau*RHS);
        %uh2=(1-Q*tau)*uh1-tau*(M\S)*uh1+tau*M\RHS;

        %uh2=(1-Q*tau)*uh1-tau*M\(S*uh1)+tau*M\RHS;
        
        %plot
        if (ti==1 || ti==3)
            lmbd=figure;
            plotfunction(uh2,T,DoFHandler); 
            title(['m=',num2str(m_size),',t=',num2str(ti),',tau=',num2str(tau),',tau size=',num2str(tau_size)]);
            %saveas(lmbd, ['C:\Users\Sunlu\Desktop\HM_610\hm4_pro1_3\hm4_prob3_ex_m',...
            %    num2str(m_size),'_t',num2str(ti),'_tau',num2str(tau_size),'.png'],'png');
            %if (ti==3 & tau_size==82)
            %    uh2
            %end
                
        end
        uh1=uh2;
    end
end
end