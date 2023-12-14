function [PM,PM_S1,PM_S2,NPM]=Model_Equations_Parasite_Mortality(t,A,B)
% Evaluates the mortality induced by parasites at any given time for the
% polymorphic population.
global m_p D m m_S m_T I N
n=length(B);

options = optimset('Display','off');
SS=fsolve(@(F)M(F,A,B),(N/n).*ones(1,n-1),options);
SS(SS<0)=0;
SS(n)=N-sum(SS);
T=SS;
S=A'-SS;

A=A';

PM=NaN(1,n);   PM_S1=PM;   PM_S2=PM;   NPM=PM;

for i=1:n
    PM(i)=(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i))+m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
    PM_S1(i)=(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i));
    PM_S2(i)=m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
end

end
