function ode=Model_Equations_Combine(t,A,B)
% Evaluates the polymorphic model as described in the main text.
global b m_p D m m_S m_T I N KG
n=length(B);

options = optimset('Display','off');
SS=fsolve(@(F)M(F,A,B),(N/n).*ones(1,n-1),options);
SS(n)=N-sum(SS);
T=SS;
S=A'-SS;

ode=NaN(n,1);
A=A';
N_T=sum(T); N_S=sum(S);

for i=1:n
    ode(i)=b(N_T+N_S)*sum(KG(i,:).*(A))-(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i))-m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
end
end
