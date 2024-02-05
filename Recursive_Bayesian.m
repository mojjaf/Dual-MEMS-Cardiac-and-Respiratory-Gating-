function [p_x_condition, E_x, Var_x]=Recursive_Bayesian(x,p_x)
%Input
%x: range of x equally spaced vector of size (1*N)
%p_x: probability density matrix of size (n*N) where p_x(i,:) is the probability density function of the ith observation
%Output
%p_x_condition: Conditional probability density function of x of size (1*N)
%E_x: expected value of x
%Var_x: variance of x
dx=x(2)-x(1);                                                     %dx used for integration
n=size(p_x,1);
p_x_condition=p_x(1,:);                                           %Initial
for i=2:1:n
    
p_x_condition_i_1=p_x_condition;                                  %From previous step
p_condition=trapz(p_x(i,:).*p_x_condition_i_1)*dx;                %p(x_i=x_i-1,.......x_2=x_1)
p_x_condition=p_x(i,:).*p_x_condition_i_1/p_condition;            %p(x=x_i|x_i=x_i-1,....x_2=x_1)=p(x=xi)*p(x=x_i-1|x_i-1=x_i-2,....x_2=x_1)/p(x_i=x_i-1,.......x_2=x_1) [Bayes Theorem]
end
E_x=trapz(x.*p_x_condition)*dx;                                   %E(x=x_n|x_n=x_n-1,....x_2=x_1) [Conditional Expectation]
Var_x=trapz((x-E_x).^2.*p_x_condition)*dx;                        %Variance of x [Conditional]     
end