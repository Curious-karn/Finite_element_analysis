
function [result]=gauss_quad(f,oder,lower_lim,upper_lim)

% function for gauss quadrature integral
 N=oder;
 a=lower_lim;
 b=upper_lim;
 [nodes,weight]=g_int(N,-1,1);
 
 x1=-a+b;
 x2=a+b;
 [nodes,weight]=g_int(N,-1,1);
 w=weight;
 psi=nodes;
 A=0;
 for i=1:length(psi)
    w_phys(i)=(a*0.5*(1-psi(i)))+(b*0.5*(1+psi(i)));
    A=A+subs(f,w_phys(i))*(w(i));

 end
 
 
 B=A*(b-a)*0.5;
 result=B;
end