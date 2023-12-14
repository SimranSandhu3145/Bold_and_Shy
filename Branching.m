function B=Branching(Nn)
%Evidence of branching behaviour in monomorphic population.
global KG threshold T

figure(3);  colormap(flipud(gray));   set(gca,'ColorScale','log')

ad_end=5000;

threshold=1e0;

all=linspace(0,1,100);
B=zeros(length(all),ad_end);
B_in=0.4;
in=find(abs(all-(B_in))==(min(abs(all-(B_in)))));
in_r=in-2:in;
Parameters('N',Nn,in_r)
KG=diag(ones(1,length(in_r)));
[~,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,all(in_r)),[0 100],(Nn+1).*ones(1,length(in_r)));
B(in_r,2)=sol(end,:);
in_r=in_r(sol(end,:)>threshold);

for i=2:ad_end
    %% Generate new mutants
in_m=min(length(all),max(0,in_r+round(20.*(2.*(rand(1,length(in_r))-0.5))))); 
in_sim=unique([in_r,in_m]);
while length(in_sim)<=length(in_r)
    in_m=min(length(all),max(0,in_r+round(20.*(2.*(rand(1,length(in_r))-0.5)))));
    in_sim=unique([in_r,in_m]);
end
%% Simulate their behaviours
if i>2
    if length(in_r)<2
        KG=diag(ones(1,length(in_sim))); %Clonal reproduction
        [t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,all(in_sim)),[0 100],(Nn+1).*ones(1,length(in_sim)));
    else
    KG=diag(ones(1,length(in_sim)));
    [t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,all(in_sim)),[0 100],B(in_sim,i)+threshold);
    end
else
    KG=diag(ones(1,length(in_sim)));
[t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,all(in_sim)),[0 100],[max(sol(end,:)),threshold.*ones(1,length(in_sim)-1)]);
end

%% Keep only the surviving strains
in_r=(in_sim(sol(end,:)>threshold));
KG=diag(ones(1,length(in_r)));
if length(in_r)>1
    T=(Nn/length(in_r)).*ones(1,length(in_r));
[t,sol]=ode45(@(t,F)Model_Equations_Combine(t,F,all(in_r)),[0 10000],sol(end,sol(end,:)>threshold));
else
    [t,sol]=ode45(@(t,F)Model_Equations_One_Strain(t,F,all(in_r)),[0 100],sol(end,sol(end,:)>threshold));
end
for j=1:length(in_r)
B(in_r,i+1)=sol(end,:);
end
if sum(sol(end,:)>threshold)<6
    in_r=(in_r(sol(end,:)>threshold))
else
    keep=nonzeros((sol(end,:)>threshold).*[1:length(sol(end,:))]);
    in_r=unique([in_r(sol(end,:)==max(sol(end,keep(1:2)))),in_r(sol(end,:)==max(sol(end,keep(3:end-1)))),in_r(sol(end,:)==max(sol(end,nonzeros(keep(3:end-1)'.*(sol(end,keep(3:end-1))<max(sol(end,keep(3:end-1)))))))),in_r(keep(end))])
end

%% Plot results
if mod(i,100)==0
    figure(3)
    contourf(all,1:i,B(:,1:i)',400)
    ylabel('Evolutionary Time');    xlabel('B');
    set(gca,'FontSize',30)
end
end

end