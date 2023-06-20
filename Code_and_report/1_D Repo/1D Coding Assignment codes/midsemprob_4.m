%% Midesem Code 4 
%  Problem No 4  ( Solve Given problem -- A cantilever subjected to a
%  continuous axial body force f(x) and a given point load (P1=5) at X=0.5L and
%  another point load at (x=L) Pl=10
%  U0=0;
% L=1; Ea=1; u0=0; x0=0; xl=L; k=0; f=10; Pl=10; P0=0; P1=5; f=10; N=[3 5 9];

clear all


%% Initiating variables

clc
x=sym('x','real');
L=1;
Ea=1;
u0=0;
x0=0;
xl=L;
k=0;
f=10;
Pl=10;
P0=0;
P1=5;
f=10;
N=[3 5 9];
X=x0:0.01*L:L;



%% Exact Solution
X1_=x0:0.01*L:(1-0.01);
syms U_e1(x)

ode1=diff(U_e1,x,2)==-f/Ea;
DU_e1=Ea*diff(U_e1,x);
cond1 = U_e1(0) == 0;
cond2 = DU_e1(0) == (Pl+P1+int(f,x,0,L));
% Exact solution for U_ex(1)
conds1=[cond1 cond2];
U_eSol1=dsolve(ode1,conds1);
U_e1=simplify(U_eSol1);
U_e1_1=double(subs(U_e1,x,0.5*L));

syms U_e2(x)

ode1=diff(U_e2,x,2)==-f/Ea;
DU_e2=Ea*diff(U_e2,x);
cond1 = U_e2(0.5*L) == U_e1_1;
cond2 = DU_e2(L) == Pl;
% Exact solution for U_ex(1)
conds2=[cond1 cond2];
U_eSol2=dsolve(ode1,conds2);
U_e2=simplify(U_eSol2);

% Making piecewise Uexact 
U_ep=piecewise(((x>=0)&(x<=0.5*L)),U_e1,((x>0.5*L)&(x<=L)),U_e2);
DU_ep=diff(U_ep);
U_e=double(subs(U_ep,x,X));
DU_e=double(subs(DU_ep,x,X));


%% All problems solution for different values of N

for i_gc=1:length(N)
  

    %% code for polynomial phi
    for i=i_gc                               %Loop for various values of N
    phi=sym(zeros(N(i),1));
    K=zeros(N(i));                                  % Empty K and F matrix
    F=zeros(N(i),1);

    % defining function for phi 
    for j=1:N(i)                                    % Loop for defining function phi
        phi(j)=x^(j-1);

    % Defining K and F matrix
        for l=1:N(i)
            F(l)=int(f*phi(l),x,0,L)+(Pl*subs(phi(l),x,L))-(P0*subs(phi(l),x,0))+(P1*subs(phi(l),x,0.5*L));
       
            for h=1:N(i)
                K(l,h)=int((Ea*diff(phi(l))*diff(phi(h))),x,0,L);

            end
        end
    end
    
    %% Tempering for end conditions
    K(1,1)=1;
    F(1)=u0;
    for a=2:N(i)
        K(1,a)=0;
        F(a)=F(a)-u0*K(a,1);
        K(a,1)=0;
    end
    
    % solution of matrix for first type of basis function

    alpha=K\F;
    U_ = ((alpha)')*phi;
    U(:,i,1)=double(subs(U_,x,X));

    % Force calculation for first type basis function
    DU=Ea*diff(U_,x);
    D_U(:,i,1)=double(subs(DU,x,X));

    % error calculation for first type of basis function
    e_N=U_ep-U_;
    B=int((Ea*diff(e_N)*diff(e_N)),x,0,L);
    err(1,i)=double(sqrt(B));
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %"code for polynomial Basis function"%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code for Hat basis function (bsis function used in problem 2)


    for i=1:i_gc                               %Loop for various values of N
    
    %Defining function phi
    phi=sym(zeros(N(i),1));
    for v=1:N(i)
        s=x;
    
        if ((v~=1)&&(v~=N(i)))
            g=1/(v-1);
            for m=1:(v-1)
                s=s*(((g*m)*L)-x);
            end

            phi(v)=s;
        
        elseif (v==1)
            phi(v)=(1-(x/L));
        else
            phi(v)=(x/L);
        end
    end

    alpha=zeros(N(i));
    K=zeros(N(i));                              % Empty K and F matrix
    F=zeros(N(i),1);

    % Defining K and F matrix
        for l=1:N(i)
            F(l)=int(f*phi(l),x,0,L)+(Pl*subs(phi(l),x,L))-(P0*subs(phi(l),x,0))+(P1*subs(phi(l),x,0.5*L));
       
            for h=1:N(i)
                K(l,h)=int((Ea*diff(phi(l))*diff(phi(h))),x,0,L);

            end
        end
       


    %% Tampering stiffness matrix 
    % for end condition (u=0 at x=0) and (u=0 at x=L) 
    % step-1:putting K(1,1)=1 and other elements of row 1 = 0 and tempering
    % F vector acordingly
    K(1,1)=1;
    F(1)=u0;
    for a=2:N(i)
        K(1,a)=0;
        F(a)=F(a)-u0*K(a,1);
        K(a,1)=0;
    end





    %% solution of matrix for first type of basis function

    alpha=K\F;
    U_ = ((alpha)')*phi;
    U(:,i,2)=double(subs(U_,x,X));

    % Force calculation for first type basis function
    DU=Ea*diff(U_,x);
    D_U(:,i,2)=double(subs(DU,x,X));

    % error calculation for first type of basis function
    e_N=U_ep-U_;
    B=int((Ea*diff(e_N)*diff(e_N)),x,0,L);
    err(2,i)=double(sqrt(B));



    end

%%%%%%%%%%%-----------------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp("Code for 2nd type of basis function ends here");
%%%%%%%%%%%------------------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %% Array of values of U at different positions of x
%     for s=1:length(X)
%         Uex(s)=double(subs(U_e,x,X(s)));
%         DUex(s)=double(subs(DU_e,x,X(s)));
%     end

    for z=i_gc
        %% Discretisation of Domain
        X2=x0:L/N(z):xl;
        for n=1:(length(X2)-1)
            elem(n,:)=[X2(n),X2(n+1)];
        end
        N_1=length(X2);    % no of nodes

        %% Weak form
        %fl= local shape function
        %fg= global shape function
        F=zeros(N_1,1);               % global load vector
        K=zeros(N_1,N_1);               % Global stiffnes matrix
        K_l=zeros(2,2,(N_1-1));            % local stiffness matrix
        F_l=zeros(2,(N_1-1));              % local load vector

        for a=1:N(z)
            L_s(a,1)=(x-elem(a,2))/(elem(a,1)-elem(a,2));   %local shape function
            L_s(a,2)=(x-elem(a,1))/(elem(a,2)-elem(a,1));   % local shape function
            DL_s=diff(L_s);
            for i=1:2
                F_l(i,a)=F_l(i,a)+int((f*L_s(a,i)),x,elem(a,1),elem(a,2));
                ig_ind=a+i-1;                                % index generation for global matrix
                F(ig_ind,1)=F(ig_ind,1)+F_l(i,a);
                for j=1:2
                    K_l(i,j,a)=K_l(i,j,a)+int(((Ea*DL_s(a,i)*DL_s(a,j))),x,elem(a,1),elem(a,2));
                    jg_ind=a+(j-1);                          % index genration for global matrix
                    K(ig_ind,jg_ind)=K(ig_ind,jg_ind)+K_l(i,j,a);
                end
            end

        end
       

  %%
  K(1,1)=1;
  F(1)=u0;
  for b=2:N_1
      K(1,b)=0;
      F(b)=F(b)-u0*K(b,1);
      K(b,1)=0;
  end
  F(1)=F(1)-P0;
  F(N_1)=F(N_1)+Pl;
  smr=((N_1)/2);
  F(smr)=F(smr)+P1;

        

   %% solution of matrix for first type of basis function

    alpha=K\F;
    %% piecewise displacement function
    U_n=0;
    for g=1:(N_1-1)
        alpha_l(g,:)=[alpha(g),alpha(g+1)];
        U_loc(g)=alpha_l(g,1)*L_s(g,1)+alpha_l(g,2)*L_s(g,2);
        U_n=piecewise((elem(g,1)<=x)&(x<=elem(g,2)),U_loc(g),U_n);
    end
    w(:,1)=subs(U_n,x,X);
    U(:,z,3)=w(:,1) ;  % array of exact solution
    DU=Ea*diff(U_n);
    D_U(:,z,3)= subs(DU,x,X);
    
    %% error calculation
    err1=U_ep-U_n;
    Derr=diff(err1);
    B=int((Ea*Derr*Derr),x,0,1);
    err(3,z)=sqrt(double(B));


    
    



    end


end

% Ploting Graph

%% Displacement plot

plot(X,U(:,1,1),'g-.',X,U(:,1,2),'r-',X,U(:,1,3),'c--')
hold on
plot(X,U_e,'-*' )
legend({"polynomial","hierarchy","Hat","Uexact"},'Location','northwest');
xlabel("Position of x");
ylabel("Displacement");
title(" U vz X for N=3",'LineWidth',15);
figure
plot(X,U(:,2,1),'g-.',X,U(:,2,2),'r-',X,U(:,2,3),'c--')
hold on
plot(X,U_e,'+:')
legend({"polynomial","hierarchy","Hat","Uexact"},'Location','northwest');
xlabel("Position of x");
ylabel("Displacement");
title("U vx X for N=5",'LineWidth',15);
figure
plot(X,U(:,3,1),'g-.',X,U(:,3,2),'r-',X,U(:,3,3),'c--')
hold on
plot(X,U_e,'p-')
legend({"polynomial","hierarchy","Hat","Uexact"},'Location','northwest');
xlabel("Position of x");
ylabel("Displacement");
title(" U vs X for N=9",'LineWidth',15);

%% Error plot
figure
plot(N,log(err(1,:)),'p--',N,log(err(2,:)),'*-',N,log(err(3,:)),'o-');
legend('polynomial','hierachy','Hat');
xlabel("No. of element");
ylabel("log(||e||)");
title("Error plot");


%% Force plot
figure
plot(X,D_U(:,1,1),'g',X,D_U(:,1,2),'b--o',X,D_U(:,1,3))
hold on
plot(X,DU_e,'r-+')
legend({"polynomial","hierarchy","Hat","F exact"},'Location','northwest');
xlabel("Position of x");
ylabel("Force");
title("Force vs X Graph for N=3",'LineWidth',15);
figure
plot(X,D_U(:,2,1),'g',X,D_U(:,2,2),'b--o',X,D_U(:,2,3))
hold on
plot(X,DU_e)
legend({"polynomial","hierarchy","Hat","F exact"},'Location','northwest');
xlabel("Position of x");
ylabel("Force");
title("Force vs X Graph for N=5",'LineWidth',15);
figure
plot(X,D_U(:,3,1),'g',X,D_U(:,3,2),'b--o',X,D_U(:,3,3))
hold on
plot(X,DU_e,'--')
legend({"polynomial","hierarchy","Hat","F exact"},'Location','northwest');
xlabel("Position of x");
ylabel("Force");
title("Force Vs X Graph for N=9",'LineWidth',15);

