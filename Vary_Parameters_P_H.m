% Varies the levels of predation and paratism for the monomorphic population.
global b m_S m_T m_p I D m N threshold
threshold=1;
figure(1); clf;
range=[0 0.25 1 7];
B=0.5;
nu=1e-3;
parameter_name='N';
parameter_range=unique([0 50 100 500 1000 3000 logspace(0,log10(3000),75)]);
m=length(parameter_range);
Total_Fish=NaN(m,length(range),2);
for HP=1:2
    if HP==1
        P=1;
        p_var='H';
    else
        H=1;
        p_var='P';
    end
    for j=1:length(range)
        for i=1:length(parameter_range)
            localfcn(p_var,range(j))
            ParametersPH(parameter_name,parameter_range(i),B,P,H)
            ode_1=@(t,A)b(sum(A))*(A)-(m.*A+m_S(B).*(A-N)+m_T(B).*N)-m_p(B).*(sum(I(B,B,N).*(A-N)+D.*I(B,B,N).*(A-N)));
            [t,sol]=ode45(@(t,F)ode_1(t,F),[0 10000],parameter_range(i)+1);
            Total_Fish(i,j,HP)=sum(sol(end,sol(end,:)>threshold));
        end
    end
end
figure(1)
subplot(1,2,1)
plot(parameter_range,Total_Fish(:,:,1),'LineWidth',2);
xlabel('N (number of shelters)');    ylabel('Population Density, F^*');
legendStrings = "H=" + string(range);
legend(legendStrings)
set(gca,'FontSize',30)
set(gca, 'XScale', 'log')
xlim([0 3000])
subplot(1,2,2)
plot(parameter_range,Total_Fish(:,:,2),'LineWidth',2);
xlabel('N (number of shelters)');    ylabel('Population Density, F^*');
legendStrings = "P=" + string(range);
legend(legendStrings)
set(gca,'FontSize',30)
set(gca, 'XScale', 'log')
xlim([0 3000])