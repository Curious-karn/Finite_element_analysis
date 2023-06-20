%% Initiating Variables 
clc
clear all    %for Precision
x=sym('x');
L=1;         %length of bar
P0=5;        % Force at end x=l
Pl=10;       % Force at end x=0
Ea=1;        % Area of Bar and multiple with Youngs modulus
k=10;        % Stifness of rubber casing 
Nel=[2 4 8];       % Number of element 
f=10*x*x;    % distributed Force on Bar
u0=0;        % Displacement at  x==0 
ul=0;        % Displacement at  x==L
x0=0; 
xl=L;        % Domain of analysis
%flag=input("enter the problem which you want to solve \n Probelm 1: represents Bar with end loads \n Probelm 2:Repesents Bar with fixed supports ");
flag=2;

%% Switching from one problem to other 
if flag==1
    %k=input('please enter value of rubber casing and 0 for no Rubber casing ');
   % P0=input('please enter the values for P0 if you want to change');
    %Pl=input('please enter the values for Pl if you want to change');
else
    flag=2;

end

%% Exact solution for all possible case 
X=linspace(x0,xl,100);
if (flag==1)&&(k==0)
    disp("exact solution does not exist");
elseif (flag==1)&&(k~=0)
    syms U_e(x)
    DU_e=diff(U_e);
    ode=diff(U_e,x,2)==k*U_e-f;
                                               % Boundary conditions
    cond1 = DU_e(0) == P0;
    cond2 = DU_e(L) == Pl;                                  
    conds=[cond1 cond2];        
                                               % Exact solution for 
    U_eSol=dsolve(ode,conds);
    U_e=simplify(U_eSol);
    DU_e=diff(U_e);

else

    syms U_e(x)
    DU_e=diff(U_e);
    ode=diff(U_e,x,2)==(1/Ea)*((k*U_e)-f);
                                                % Boundary conditions 
    cond1 = U_e(0) == 0;
    cond2 = U_e(L) == 0;
    conds=[cond1 cond2];

                                                % Exact solution for 
    U_eSol=dsolve(ode,conds);
    U_e=simplify(U_eSol);

    DU_e=diff(U_e);




end

%% Array of values of U at different positions of x
for s=1:length(X)
    Uex(s)=double(subs(U_e,x,X(s)));
    DUex(s)=double(subs(DU_e,x,X(s)));
end

for z=1:length(Nel)
    %% Discretisation of Domain
    X1=x0:L/Nel(z):xl;
    for n=1:(length(X1)-1)
        elem(n,:)=[X1(n),X1(n+1)];
    end
    N=length(X1);    % no of nodes

    %% Weak form
    %fl= local shape function
    %fg= global shape function
    F=zeros(N,1);               % global load vector
    K=zeros(N,N);               % Global stiffnes matrix
    K_l=zeros(2,2,(N-1));            % local stiffness matrix
    F_l=zeros(2,(N-1));              % local load vector

    for a=1:Nel(z)
        L_s(a,1)=(x-elem(a,2))/(elem(a,1)-elem(a,2));   %local shape function
        L_s(a,2)=(x-elem(a,1))/(elem(a,2)-elem(a,1));   % local shape function
        DL_s=diff(L_s);
        for i=1:2
            F_l(i,a)=F_l(i,a)+int((f*L_s(a,i)),x,elem(a,1),elem(a,2));
            ig_ind=a+i-1;                                % index generation for global matrix
            F(ig_ind,1)=F(ig_ind,1)+F_l(i,a);
            for j=1:2
                K_l(i,j,a)=K_l(i,j,a)+int(((Ea*DL_s(a,i)*DL_s(a,j))+k*L_s(a,i)*L_s(a,j)),x,elem(a,1),elem(a,2));
                jg_ind=a+(j-1);                          % index genration for global matrix
                K(ig_ind,jg_ind)=K(ig_ind,jg_ind)+K_l(i,j,a);
            end
        end

    end
    %% optimization for different boundary condition
    if flag==2

        % for end condition (u=0 at x=0) and (u=0 at x=L)
        % step-1:putting K(1,1)=1 and other elements of row 1 = 0 and tempering
        % F vector acordingly
        K(1,1)=1;
        F(1)=u0;
        for b=2:N
            K(1,b)=0;
            F(b)=F(b)-u0*K(b,1);
            K(b,1)=0;
        end

        % Step-2: Similar proceedure for K(n,n)=1 and so on
        K(N,N)=1;
        F(N)=ul;
        for c=(N-1):-1:2
            K(N,c)=0;
            F(c)=F(c)-ul*K(c,N);
            K(c,N)=0;
        end
    elseif ((flag==1)&&(k==0))
        %K matrix is singular, no unique solution exists\n'
        %change Pl such that F(N)=0 and change K(N,N)=1\n
        K(:,N)=0;
        K(N,:)=0;
        K(N,N)=1;
        Pl=-F(N);
        F(1)=F(1)-P0;
        F(N)=F(N)+Pl;
    else
        F(1)=F(1)-P0;
        F(N)=F(N)+Pl;


    end

    %% Solving equation
    alpha=K\F;
    t=linspace(0,1,50);

    %% piecewise displacement function
    U=0;
    for g=1:(N-1)
        alpha_l(g,:)=[alpha(g),alpha(g+1)];
        U_loc(g)=alpha_l(g,1)*L_s(g,1)+alpha_l(g,2)*L_s(g,2);
        U=piecewise((elem(g,1)<=x)&(x<=elem(g,2)),U_loc(g),U);
    end
    w(1,:)=subs(U,x,t);
    Val(z,:)=w(1,:) ;  % array of exact solution
    DVal=Ea*diff(U);
    dw(1,:)= subs(DVal,x,t);
    Val1(z,:)=dw(1,:);
    
    %% error calculation
    err(z)=U_e-U;
    Derr=diff(err(z));
    B(z)=int(((Ea*Derr*Derr)+k*Derr*Derr),x,0,1);
    errs=B(z);
    err(z)=sqrt(double(errs));
end



%% ploting Uexact and U approx
for z1=1:length(Nel)
    
    plot(t,Val(z1,:));
    hold on
end
Val_e=subs(U_e,x,t);
plot(t,Val_e)
xlabel("Position");
ylabel("Displacement");
legend("U_app N=2","U_app N=4","U_app N=8","U exact");

%% ploting D(Uexat) and D(Uapprox)
figure
for z2=1:length(Nel)
% DVal=Ea*diff(U);
% DVal_e=Ea*diff(U_e);
% Val1= subs(DVal,x,t);
% Val_e1=subs(DVal_e,x,t);
    plot(t,Val1(z2,:));
    hold on
end
DVal_e=Ea*diff(U_e);
Val_e1=subs(DVal_e,x,t);
plot(t,Val_e1);
xlabel("Position");
ylabel("Displacement");
legend("DU_app N=2","DU_app N=4","DU_app N=8","DU exact");
title("Force Plot for");



%% error ploting
figure
plot(Nel,err);
xlabel("No of element");
ylabel("log ||error||");
title("Error plot",LineStyle="-");







