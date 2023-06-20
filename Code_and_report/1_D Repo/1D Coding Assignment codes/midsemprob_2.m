%% Initiating variable 
clc
digits(64)
clear all
Ea=1;
k=[0 10];
syms x real
f=[10,10*x,10*x*x];
L=1;
N=[2 4 6 8];
U = sym(zeros(length(N),1));
z=1;
f=10*x*x;
u0=0;
ul=0;


syms U_e(x)
DU_e=diff(U_e);
ode=diff(U_e,x,2)==(1/Ea)*((k(z)*U_e)-f);

cond1 = U_e(0) == 0;
cond2 = U_e(L) == 0;

% Exact solution for U

conds=[cond1 cond2];
U_eSol=dsolve(ode,conds);
U_e=simplify(U_eSol);
% end

%%
for i=1:length(N)                               %Loop for various values of N
    
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
            F(l)=int(f*phi(l),x,0,L); %(Pl*subs(phi(l),x,L))-P0*subs(phi(l),x,0);
       
            for h=1:N(i)
                K(l,h)=int((Ea*diff(phi(l))*diff(phi(h))+(k(z)*phi(l)*phi(h))),x,0,L);

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
  
    % Step-2: Similar proceedure for K(n,n)=1 and so on
    K(N(i),N(i))=1;
    F(N(i))=ul;
    for b=(N(i)-1):-1:2
        K(N(i),b)=0;
        F(b)=F(b)-ul*K(b,N(i));
        K(b,(N(i)))=0;
    end





    %% checking for condition of k
                        
    alpha=K^(-1)*F;
    U(i,1) = (alpha')*phi;

    %% Ploting displacement 
    t=linspace(0,L,50);
    U_ = U(i,1);
    DU_=diff(U_);
    Val=subs(U_,x,t);
    DVal(:,i)=Ea*subs(DU_,x,t);
    plot(t,Val)
    b(i)=num2str(N(i));
    
    hold on
    %% defining error in the approximate solution
   
    e_N(i)=U_e-U_;
    B(i)=int(((Ea*diff(e_N(i))*diff(e_N(i)))+(k(z)*e_N(i)*e_N(i))),x,0,L);
    err(i)=double(sqrt(B(i)));
 
    
    


end

xlabel("Position (X)");
ylabel('Displacement');
title('U vs X','FontSize',20,'Color','r');

%% Ploting exact solution 
UeVal=subs(U_e,x,t);
plot(t,UeVal,'--','linewidth',2);
legend('N=2','N=4','N=6','N=8','exact');

  


%% Ploting error
figure
plot(N,log(err));        % logerthims of error
legend('Error plot');
xlabel('Number of Iteration');
ylabel('log(||e||)');





figure
for q=1:length(N)
    
    plot(t,DVal(:,q));
    
    hold on
end
xlabel("Position (X)");
ylabel('Force ');
title('Force vs Position','FontSize',20,'Color','g');

%% Ploting exact force
DU_e=diff(U_e);
UeVal=subs(DU_e,x,t);
plot(t,UeVal,'--','linewidth',2);
legend('N=2','N=4','N=6','N=8','exact');






    

