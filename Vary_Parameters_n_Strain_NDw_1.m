% Varies both N and D_w for the polymorphic population.
global KG threshold T

threshold=1;

figure(1); clf; xlabel('N (number of shelters)');    ylabel('Population Density, F^*');    set(gca,'FontSize',30)
xlim([0 3000])
figure(2); clf;
figure1 = figure(2);
% Create textbox
annotation(figure1,'textbox',...
    [0.0648191465158607 0.888049592220557 0.02634375 0.060377358490566],...
    'String',{'(A)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.344506646515861 0.885893257988752 0.0263437499999999 0.060377358490566],...
    'String','(B)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.626798313182527 0.888049592220558 0.02634375 0.060377358490566],...
    'String','(C)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.501277479849194 0.425515899498187 0.0263437500000001 0.0603773584905662],...
    'String','(E)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.0585691465158606 0.419585980360721 0.02634375 0.0603773584905662],...
    'String','(D)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',40,...
    'FitBoxToText','off');
subplot(2,3,1)
xlabel('B');    ylabel('F/sum(F)'); title('Clonal');    set(gca,'FontSize',20); box on
ylim([0 inf])
subplot(2,3,2)
xlabel('B');    ylabel('F/sum(F)'); title('D_W=0.02');  set(gca,'FontSize',20); box on
ylim([0 inf])
subplot(2,3,3)
xlabel('B');    ylabel('F/sum(F)'); title('Uniform');   set(gca,'FontSize',20); box on
ylim([0 inf])
subplot(2,2,3); xlabel('N (number of shelters)');    ylabel({'% Sheltered Bold', 'in F^*'});   set(gca,'FontSize',20)
xlim([0 3000])
subplot(2,2,4); xlabel('N (number of shelters)');    ylabel({'% Sheltered Bold', 'in all Shelters'});   set(gca,'FontSize',20)
xlim([0 3000])
figure(3); clf;
subplot(2,2,1)
xlabel('N (number of shelters)');    ylabel('Parasite Induced Mortality');    set(gca,'FontSize',15)
xlim([0 3000])
subplot(2,2,3)
xlabel('N (number of shelters)');    ylabel({'Non-Fighting','Parasite Induced Mortality'});    set(gca,'FontSize',15)
xlim([0 3000])
subplot(2,2,4)
xlabel('N (number of shelters)');    ylabel({'Fighting','Parasite Induced Mortality'});    set(gca,'FontSize',15)
xlim([0 3000])
subplot(2,2,2)
xlabel('N (number of shelters)');    ylabel({'% Parasite Induced Mortality','caused by Fighting'});    set(gca,'FontSize',15)
xlim([0 3000])

Dw_all=[0.001 0.02 0.1 1000];
nu=0.05;
B=linspace(0,1,20);    n=length(B);
parameter_name='N';
parameter_range=unique([0 50 100 500 1000 3000 logspace(0,log10(3000),30)]);
m=length(parameter_range);
All_Fish=NaN(n,m,length(Dw_all));   All_T=NaN(n,m,length(Dw_all));   All_S=NaN(n,m,length(Dw_all));
Parasite=NaN(n,m,length(Dw_all));    Parasite_NF=NaN(n,m,length(Dw_all));    Parasite_F=NaN(n,m,length(Dw_all));
for j=1:length(Dw_all)
    Total_Fish=NaN(1,length(parameter_range));
    percent=Total_Fish;
    percent_T=Total_Fish;
    Parasite_percent=Total_Fish;
    Dw=Dw_all(j);
    if Dw<0.01
        KG=diag(ones(1,length(B)));
    else
    KG=NaN(length(B),length(B));
    for k=1:length(B)
        KG(k,:)=exp(-((B-B(k)).^2)./Dw)./sum(exp(-((B-B(k)).^2)./Dw));
    end
    end
    for i=1:length(parameter_range)
        T=(parameter_range(i)/n).*ones(1,n);
        Parameters(parameter_name,parameter_range(i),B)
        if Dw<0.01
            [t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,B),[0 30000],[ones(1,n-1),(parameter_range(i)+1)]);
        else
            [t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,B),[0 1000],(parameter_range(i)+1).*ones(1,n));
        end        
        figure(5)
        plot(t,sol)
        Total_Fish(i)=sum(sol(end,:))
        All_Fish(:,i,j)=sol(end,:)';
        if ismember(parameter_range(i),[50 100 500 1000 3000])==1
            if Dw<0.01 || Dw == 0.02 || Dw>1
                figure(2)
                if Dw<0.01
                    subplot(2,3,1)
                elseif Dw==0.02
                    subplot(2,3,2)
                else
                    subplot(2,3,3)
                end
                hold all
                plot(B,sol(end,:)./sum(sol(end,:)),'LineWidth',2)
                if Dw<0.01
                legendStrings = "N = " + string([50 100 500 1000 3000]');
                legend(legendStrings)
                end
                drawnow
            end
        end
        SS=mean(sol(end-100:end,:))'
        [MP,MP_S1,MP_S2]=Model_Equations_Parasite_Mortality(t(end),sol(end,:)',B);
        Parasite(:,i,j)=MP;
        Parasite_NF(:,i,j)=MP_S1;
        Parasite_F(:,i,j)=MP_S2;
        Parasite_percent(i)=100.*sum(MP_S2)./sum(MP);
%         
%         PC_Parasite(i,j)=sum(MP)./Total_Fish(i);
%         PC_Parasite_NF(i,j)=sum(MP_S1)./Total_Fish(i);
%         PC_Parasite_F(i,j)=sum(MP_S2)./Total_Fish(i);
        
        options = optimset('Display','off');
        SS=fsolve(@(F)M(F,sol(end,:)',B),T(1:end-1),options);
        SS(SS<0)=0;
        SS(n)=parameter_range(i)-sum(SS);
        All_T(:,i,j)=SS';
        All_S(:,i,j)=(sol(end,:)-SS)';
        
        percent(i)=100.*sum(SS(B>0.7))./Total_Fish(i);
        percent_T(i)=100.*sum(SS(B>0.7))./sum(SS);
        
    end
    figure(1)
    hold all
    plot(parameter_range,Total_Fish,'LineWidth',2);
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    drawnow
    set(gca, 'XScale', 'log')
    
    figure(2)
    subplot(2,2,3)
    hold all
    plot(parameter_range,percent,'LineWidth',2)
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    set(gca, 'XScale', 'log')
        subplot(2,2,4)
    hold all
    plot(parameter_range,percent_T,'LineWidth',2)
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
set(gca, 'XScale', 'log')
    
    figure(3)
    subplot(2,2,1)
    hold all
    plot(parameter_range,sum(Parasite(:,:,j)),'LineWidth',2);
%     legendStrings = "D_W = " + string(Dw_all(1:j)');
%     legend(legendStrings)
    legend({'Clonal','D_W=0.02','D_W=0.1','Uniform','Single Strain'}')
    set(gca, 'XScale', 'log')
    subplot(2,2,3)
    hold all
    plot(parameter_range,sum(Parasite_NF(:,:,j)),'LineWidth',2);
    set(gca, 'XScale', 'log')
    subplot(2,2,4)
    hold all
    plot(parameter_range,sum(Parasite_F(:,:,j)),'LineWidth',2);
    drawnow
    set(gca, 'XScale', 'log')
    subplot(2,2,2)
    hold all
    plot(parameter_range,Parasite_percent,'LineWidth',2);
    drawnow
    set(gca, 'XScale', 'log')
end