% Varies both N and D_w for the monomorphic population.
global KG threshold T

threshold=1;

figure(1); clf; xlabel('N (number of shelters)');    ylabel('Population Density, F^*');    set(gca,'FontSize',30)
figure(2); clf;
subplot(2,3,1)
xlabel('B');    ylabel('F/sum(F)'); title('Clonal');    set(gca,'FontSize',20); box on
subplot(2,3,2)
xlabel('B');    ylabel('F/sum(F)'); title('D_W=0.02');  set(gca,'FontSize',20); box on
subplot(2,3,3)
xlabel('B');    ylabel('F/sum(F)'); title('Uniform');   set(gca,'FontSize',20); box on
subplot(2,2,3); xlabel('N (number of shelters)');    ylabel('% Sheltered Bold in F^*');   set(gca,'FontSize',20)
subplot(2,2,4); xlabel('N (number of shelters)');    ylabel('% Sheltered Bold in all Shelters');   set(gca,'FontSize',20)
figure(3); clf;
subplot(2,2,1:2)
xlabel('N (number of shelters)');    ylabel('Parasite Induced Mortality');    set(gca,'FontSize',15)
subplot(2,3,4)
xlabel('N (number of shelters)');    ylabel({'Non-Fighting','Parasite Induced Mortality'});    set(gca,'FontSize',15)
subplot(2,3,5)
xlabel('N (number of shelters)');    ylabel({'Fighting','Parasite Induced Mortality'});    set(gca,'FontSize',15)
subplot(2,3,6)
xlabel('N (number of shelters)');    ylabel({'% Parasite Induced Mortality','caused by Fighting'});    set(gca,'FontSize',15)

nu=0.05;
B_all=[0.5];    n=1;
parameter_name='N';
parameter_range=unique([0 50 100 500 1000 3000 logspace(0,log10(3000),30)]);
m=length(parameter_range);
All_Fish=NaN(n,m,1);   All_T=NaN(n,m,1);   All_S=NaN(n,m,1);
Parasite=NaN(n,m,1);    Parasite_NF=NaN(n,m,1);    Parasite_F=NaN(n,m,1);
for j=1:length(B_all)
    B=B_all(j)
    Total_Fish=NaN(1,length(parameter_range));
    percent=Total_Fish;
    percent_T=Total_Fish;
    Parasite_percent=Total_Fish;
    for i=1:length(parameter_range)
        Parameters(parameter_name,parameter_range(i),B)
        [t,sol]=ode45(@(t,F)Model_Equations_One_Strain(t,F,B),[0 500],(parameter_range(i)+1).*ones(1,n));
        Total_Fish(i)=sum(sol(end,:))
        
        [MP,MP_S1,MP_S2]=Model_Equations_Parasite_Mortality_one_strain(t(end),sol(end,:)',B);
        Parasite(:,i,j)=MP;
        Parasite_NF(:,i,j)=MP_S1;
        Parasite_F(:,i,j)=MP_S2;
        Parasite_percent(i)=100.*sum(MP_S2)./sum(MP);
        
        
        percent(i)=100.*sum(0)./Total_Fish(i);
        percent_T(i)=100.*sum(0)./sum(1);
        
    end
    figure(1)
    hold all
    plot(parameter_range,Total_Fish,'--','LineWidth',2);
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    drawnow
    set(gca, 'XScale', 'log')
    
    figure(2)
    subplot(2,2,3)
    hold all
    plot(parameter_range,percent,'--','LineWidth',2)
%     plot(parameter_range,100.*(Total_Fish-parameter_range)./Total_Fish,'LineWidth',2)
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    set(gca, 'XScale', 'log')
        subplot(2,2,4)
    hold all
    plot(parameter_range,percent_T,'--','LineWidth',2)
%     plot(parameter_range,100.*(Total_Fish-parameter_range)./Total_Fish,'LineWidth',2)
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
set(gca, 'XScale', 'log')
    
    figure(3)
    subplot(2,2,1:2)
    hold all
    plot(parameter_range,(Parasite(:,:,j)),'--','LineWidth',2);
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    set(gca, 'XScale', 'log')
    subplot(2,3,4)
    hold all
    plot(parameter_range,(Parasite_NF(:,:,j)),'--','LineWidth',2);
    set(gca, 'XScale', 'log')
    subplot(2,3,5)
    hold all
    plot(parameter_range,(Parasite_F(:,:,j)),'--','LineWidth',2);
    drawnow
    set(gca, 'XScale', 'log')
    subplot(2,3,6)
    hold all
    plot(parameter_range,Parasite_percent,'--','LineWidth',2);
    drawnow
    set(gca, 'XScale', 'log')
end