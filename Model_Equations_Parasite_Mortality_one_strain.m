function [PM,PM_S1,PM_S2,NPM]=Model_Equations_Parasite_Mortality_one_strain(t,A,B)
% Evaluates the mortality induced by parasites at any given time for the
% monomorphic population.
global m_p D m m_S m_T I N
n=length(B);

T=N;
S=A'-T;

PM=NaN(1,n);   PM_S1=PM;   PM_S2=PM;   NPM=PM;
A=A';

for i=1:n
    PM(i)=(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i))+m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
    PM_S1(i)=(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i));
    PM_S2(i)=m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
end

end
