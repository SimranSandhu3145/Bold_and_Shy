function ode=Model_Equations_One_Strain(t,A,B)
% Evaluates the monomorphic population as described in the main text.
global b m_p D m m_S m_T I N
n=length(B);

T=N;
S=A'-T;

ode=NaN(n,1);
A=A';
N_T=sum(T); N_S=sum(S);

for i=1:n
    ode(i)=b(N_T+N_S)*(A)-(m.*A(i)+m_S(B(i)).*S(i)+m_T(B(i)).*T(i))-m_p(B(i)).*(sum(I(B(i),B,T).*S(i)+D.*I(B,B(i),T(i)).*S));
end
end
