function [DU_s,a_s,K_s,m] = scr(Nel,order_P_c,DUfem,DU_exact,Length_of_domain)

% function for Super-convergence rate

              % order of lagrangian basis function 
n=Nel;          % number of element
p_c=order_P_c;        % order of curve_fit
L=Length_of_domain;
p_m=p_c+1;
x0=0;
xl=L;
Ea=1;
syms x
x1=linspace(x0,xl,n+1);

for l=1:n
    m_bound(l,:)=[x1(l),x1(l+1)];         % upper and lower limits for each mesh element [lw,up]
    m(l,:)=g_int(p_c,x1(l),x1(l+1));      % Gauss_quad points for given order of curve fit

end

% initialisation for K_tilda and F_tilda matrices
K_s=zeros((p_m),(p_m));
F_s=zeros((p_m),1);
a_s=zeros((p_m),1);
DU=DUfem;
DU_e=DU_exact;

%% Gauss integration points in each element
for k=1:n
    a=m_bound(k,1);
    b=m_bound(k,2);
    g_pnt(k,:)=g_int(p_c,a,b);
end

%%  
for k=1:n
    if k==n
            l=k;          % for first boundary elements
             for i=1:p_m
                 f_s1=(DU*x^(i-1));
                 od=5;
                 F_s(i,1)=gauss_quad(f_s1,od,m_bound(l,1),m_bound(l,2))+gauss_quad(f_s1,od,m_bound((l-1),1),m_bound((l-1),2));
                 for j=1:p_m
                     od2=6;
                     k_s=(x^(i-1))*(x^(j-1));
                     K_s(i,j)=gauss_quad(k_s,od,m_bound(l,1),m_bound(l,2))+gauss_quad(k_s,od,m_bound((l-1),1),m_bound((l-1),2));
                 end
             end
             a_s(:,l)=K_s\F_s;
        
    elseif k==1
            l=k;              % for first boundary element
               for i=1:p_m
                   f_s1=DU*x^(i-1);       % gauss_quad functions for element-(L)
                   od=5;
                   F_s(i,1)=gauss_quad(f_s1,od,m_bound(l,1),m_bound(l,2))+gauss_quad(f_s1,od,m_bound((l+1),1),m_bound((l+1),2));
                 for j=1:p_m
                     od2=6;
                     k_s=(x^(i-1))*(x^(j-1));
                      K_s(i,j)=gauss_quad(k_s,od,m_bound(l,1),m_bound(l,2))+gauss_quad(k_s,od,m_bound((l+1),1),m_bound((l+1),2));
                 end
               end
               a_s(:,l)=K_s\F_s;
    else   
        for l=2:(n-1)              % for inner elements
            for i=1:p_m
                f_s1=DU*x^(i-1);
                od=5;
                F_s(i,1)=gauss_quad(f_s1,od,m_bound(l,1),m_bound(l,2))+gauss_quad(f_s1,od,m_bound((l-1),1),m_bound((l-1),2))+gauss_quad(f_s1,od,m_bound((l+1),1),m_bound((l+1),2));
                for j=1:p_m
                    od=6;
                    k_s=(x^(i-1))*(x^(j-1));
                    K_s(i,j)=gauss_quad(k_s,od,m_bound(l,1),m_bound(l,2))+gauss_quad(k_s,od,m_bound((l-1),1),m_bound((l-1),2))+gauss_quad(k_s,od,m_bound((l+1),1),m_bound((l+1),2));
                end
            end
            a_s(:,l)=K_s\F_s;
        end
    end
end
%%
DU_s=0;
for l=1:n
    for i=1:p_m
        DU_sml(l,i)=a_s(i,l)*x^(i-1);
    end
    DU_sm(l,1)=sum(DU_sml(l,:));
    DU_s=piecewise((x>=m_bound(l,1))&(x< m_bound(l,2)),DU_sm(l),DU_s);

end
%%

t=linspace(0,1,100);
val_s=subs(DU_s,t);
val=subs(DU,t);
val_e=subs(DU_e,t);
hold on
plot(t,val_s,'g-')
plot(t,val,'r-')
plot(t,val_e,'--',LineWidth=1.5);
hold off

%%
legend('SCR','FEM','Exact');


end