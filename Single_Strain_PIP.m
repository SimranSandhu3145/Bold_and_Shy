function Single_Strain_PIP(Nn)
%Construction of the PIP for the monomorphic population based on the
%adaptive dynamics framework as outlined in SM6.
figure(2); clf;
global b m_p D m m_S m_T N b0 k nu omega
Parameters('N',Nn,0.5)

B_all=linspace(0,1,500);
PIP=NaN(length(B_all));

for i=1:length(B_all)
    Br=B_all(i);
    
    sts=roots([-b0./k b0-m-m_S(Br)-(1+D)*m_p(Br)*nu(Br)*N (m_S(Br)-m_T(Br))*N+(1+D)*m_p(Br)*nu(Br)*N^2]);
    sts=min(nonzeros(sts.*(sts>0)));
    if length(sts)>1 
        pause
    end
    for j=1:length(B_all)
        Bm=B_all(j);
        Tm=N*nu(Bm)*omega(Bm,Br)/((sts-N)*nu(Br)*omega(Br,Bm)+nu(Bm)*omega(Bm,Br)*N);
        lambda=b(sts)-m-m_S(Bm)+(m_S(Bm)-m_T(Bm)-2*omega(Br,Bm)*(m_p(Bm)+m_p(Br)*D)*(sts-N)*nu(Br))*Tm;        
        
        if lambda>0
            PIP(i,j)=0;
        else
            PIP(i,j)=1;
        end
    end
    
end
figure(2)
colormap(flipud(gray));
hold all
contourf(B_all,B_all,PIP',1)
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
xlabel('B_r');    ylabel('B_m');
set(gca,'FontSize',30)
box on
drawnow