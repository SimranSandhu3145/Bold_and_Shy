function M=M(T,A,B)
% function to evalute the function M, describing the fast shelter-shoal
% transitioning actions. The function is described in detail in the main
% text.
global N omega nu

T(T<0)=0;

n=length(B);
M=NaN(1,n-1);
for i=1:n-1
    M(i)=(A(i)-T(i)).*(sum(nu(B(i)).*omega(B(i),B(1:end-1)).*T)+nu(B(i)).*omega(B(i),B(end)).*(N-sum(T)))-sum((A(1:end-1)'-T).*nu(B(1:end-1)).*omega(B(1:end-1),B(i)).*T(i))-(A(end)-(N-sum(T))).*nu(B(end)).*omega(B(end),B(i)).*T(i);
end