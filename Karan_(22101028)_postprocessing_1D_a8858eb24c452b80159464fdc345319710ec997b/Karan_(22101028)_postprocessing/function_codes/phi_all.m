function [shape_fun,loc_node,loc_inelem] = phi_all(Oder_p,num_elem_N,length_of_domain)


% this function will return the shape functions for all elements 

%% function for piecewise function 
syms x

%% Preprocessor 
L=length_of_domain;
N=num_elem_N;      % Number of elements
h=L/N;             % Uniform Mesh size 
p=Oder_p;          % order of basis function
n=p+1;             % Number of Nodes inside an element
n_dom=(N*p)+1;     % Number of nodes in Domain (0,L)
n_loc=linspace(0,1,((N*p)+1));  % Locations of nodes inside the domain
x_jun=0;

% storing elementwise location
loc_c=2;           % counter of location in domain

% with the help of below loops we are storing location of each element in a
% perticuler row

for i=1:N           % for element run
    for j=1:(p+1)   % for nodes inside an element
        if j==1
            loc(i,j)=x_jun;
        else
        loc(i,j)=n_loc(loc_c);
        loc_c=loc_c+1;
        end
        
    end
    x_jun=loc(i,(p+1));
end
 

% shape function
% up_t=1;     % initialisation of upper term 
% lw_t=1;     % initialisation of lower term
for i1=1:N                    
    for k1=1:p+1
        up_t=1;     % initialisation of upper term 
        lw_t=1;     % initialisation of lower term

        for j1=1:(p+1)
            if j1~=k1
                up_t=up_t*(x-loc(i1,j1)); % upper term of shape function
                lw_t=lw_t*(loc(i1,k1)-loc(i1,j1)); % lower term of shape function 
            end
        end
        shape_f(k1,i1)=(up_t)/(lw_t);  % final shape function
        % Note: Here shape functions are being stored in a coloumn   
        % (shape function,k_th_elem)
    end

end
% shape_ft=shape_f';
% % for storing shape function in a list
% c=1;
% for i2=1:N
%     for j2=1:p
%         shape_fl(c)=shape_f(j2,i2);
%         c=c+1;
%     end
% end
% shape_fl(c)=shape_f(p+1,N);
% 
%     
% lshape_fun=shape_fl;  
shape_fun=shape_f;      %Returning the result
loc_node=n_loc;
loc_inelem=loc;


%% ploting of shape functions
for l=1:N
    t(:,l)=linspace(((L/N)*(l-1)),((L/N)*l),100);
end
%%

for i2=1:N
    for j2=1:p+1
        f1=shape_f(j2,i2);
        val=subs(f1,x,t(:,i2));
        plot(t(:,i2),val);
        hold on
        grid on
    end
end
y=zeros(100,1);
plot(t,y,'r--');
ylim([-0.4,1.5]);
Na=strcat("N = ",num2str(N)," p = ",num2str(p));
[q,s]=title("Shape function plot for N=",Na);


end