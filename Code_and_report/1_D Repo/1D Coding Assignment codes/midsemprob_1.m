%% Initiation of variables
clc
clear all
Ea=1;
k=[0 10];
Pl=10;
P0=5;
x=sym('x','real');
f=10*x*x;
L=1;
N=[2 4 6 8];
U = sym(zeros(length(N),1));
z=1;

if k(z)==0
    disp('exact solution is not defined');

else

    syms U_e(x)
    DU_e=diff(U_e);
    ode=diff(U_e,x,2)==k(z)*U_e-f;
    cond1 = DU_e(0) == P0;
    cond2 = DU_e(L) == Pl;
    
   % Exact solution for U

    conds=[cond1 cond2];
    U_eSol=dsolve(ode,conds);
    U_e=simplify(U_eSol);
end



% vectors for storing errors
B=zeros(length(N),1);
err=zeros(length(N),1);

%% process of code


for i=1:length(N)                               %Loop for various values of N
    phi=sym(zeros(N(i),1));
    K=zeros(N(i));                         % Empty K and F matrix
    F=zeros(N(i),1);

    % defining function for phi 
    for j=1:N(i)                                % Loop for defining function phi
        phi(j)=x^(j-1);

    % Defining K and F matrix
        for l=1:N(i)
            F(l)=int(f*phi(l),x,0,L)+(Pl*subs(phi(l),x,L))-P0*subs(phi(l),x,0);
       
            for h=1:N(i)
                K(l,h)=int((Ea*diff(phi(l))*diff(phi(h))+(k(z)*phi(l)*phi(h))),x,0,L);

            end
        end
    end
    
    %% checking for condition of k
    if k(z)~=0                           % if k is non zero the this if block will execute 
        alpha=K^(-1)*F;
        U(i,1) = (alpha')*phi;


    else                              % if k is zero then execuute this block of code
        Pl=P0*subs(phi(1),x,0)-int(f*phi(1),x,0,L);
        for n=1:N(i)
            F(n)=int(f*phi(n),x,0,L)+(Pl*subs(phi(n),x,L))-P0*subs(phi(n),x,0);
            for m=1:N(i)
                K(n,m)=int((Ea*diff(phi(n))*diff(phi(m))),x,0,L);
                K(1,1)=1;
            end
        end
        alpha=K^(-1)*F;
        U(i,1) = (alpha')*phi;
    end

        
    t=linspace(0,L,50);
    U_ = U(i,1);
    DU_(i)=diff(U(i,1));
    Val=subs(U_,x,t);
    plot(t,Val)
    b(i)=num2str(N(i));
    
    hold on

%% defining error in the approximate solution
    if k(z)~=0
        e_N(i)=U_e-U_;
    
        B(i)=int(((Ea*diff(e_N(i))*diff(e_N(i)))+(k(z)*e_N(i)*e_N(i))),x,0,L);
        err(i)=log(double(sqrt(B(i))));
 
    
    end
    



end  

legend('N=2','N=4','N=6','N=8');
xlabel("Position (X)");
ylabel('Displacement');
title('U vs X','FontSize',20,'Color','r');

%% Ploting exact solution 
if k(z)~=0
    UeVal=subs(U_e,x,t);
    plot(t,UeVal,'--','linewidth',2);
    legend('N=2','N=4','N=6','N=8','exact');



    %% Ploting error
    figure
    plot(N,err);
    legend('Error plot');
    xlabel('Order of polynomial basis');
    ylabel('log ||e||');
end

%% Ploting force
figure
for m=1:length(N)
    Val2=subs(DU_(m),x,t);
    Val_2=double(Val2);
    plot(t,Val_2);
    hold on
end
% DU_e=diff(U_e);
% UeVal1=subs(U_e,x,t);
% plot(t,UeVal1);
title("Force plot");
xlabel("position");
ylabel("Force");
legend('N=2','N=4','N=6','N=8');



