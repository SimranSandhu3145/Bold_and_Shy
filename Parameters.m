function Parameters(parameter_name,parameter_range,B)
%Parameter values used 
global b m_p D m m_S m_T I N b0 k nu omega KG

N=200;
b0=0.15;
k=10000;
b=@(x)b0.*(1-x/10000);


epsilon=0.1;
m_B=@(B)1-epsilon*B;

gamma1S=0.02;
gamma2S=0.05;
gamma1T=0.035;
gamma2T=0;
gamma_H=1;
gamma_P=0.01;

H=1;    P=1;
m_S=@(B)m_B(B)*(gamma1S*H+gamma2S*P);
m_T=@(B)m_B(B)*(gamma1T*H+gamma2T*P);
m_p=@(B)m_B(B)*(gamma_H*H);
m=gamma_P*P;


D=0.1;

localfcn(parameter_name,parameter_range)


dw=30;
omega=@(B_i,B_j)(exp(-dw.*(B_j-B_i))./(1+exp(-dw.*(B_j-B_i))));
Bv=0.5;
m2=5;
nu_m=1e-3;
nu=@(B)nu_m.*B.^m2./(B.^m2+Bv.^m2);
I=@(B_i,B_j,T2)2*nu(B_i).*(omega(B_i,B_j)).*T2;
end