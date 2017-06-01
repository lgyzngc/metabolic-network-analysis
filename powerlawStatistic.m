filenameList = dir('D:\huansheng\program\biggstat\format');
network_statisticTable={};
cofactorFreeDegree=[];
GroupreactionDegree=[];
index=1;
for i=3:size(filenameList,1)
    network_statisticTable{index} = csvread(strcat('D:\huansheng\program\biggstat\format\',filenameList(i).name),1,0);
    str=regexp(filenameList(i).name,'groupreaction','match');
    if isempty(str)
        cofactorFreeDegree=[cofactorFreeDegree,index];
    else
        GroupreactionDegree=[GroupreactionDegree,index];
    end
    index=index+1;
end

dataorgin=load('powerlawdata.txt');
cluster_coefficient_cofactorfree={};
cluster_coefficient_groupreaction={};

for i=1:length(cofactorFreeDegree)
    tempdata=network_statisticTable{cofactorFreeDegree(i)};
    cluster_coefficient_cofactorfree{i}=tempdata(:,2);
    network_degreeFreq=tabulate(tempdata(:,5));
    temp=[];
    index=1;
    for j=1:size(network_degreeFreq,1)
        if network_degreeFreq(j,1)~=0
            temp(index,:)=network_degreeFreq(j,:);
            index=index+1;
        end
    end
    network_degreeFreq=temp;
    network_degreeFreq=log(network_degreeFreq(:,[1,2]));
    degreecutoff=find(network_degreeFreq(:,1)<=3.4);
    network_degreeFreq=network_degreeFreq(degreecutoff,:);
    if i==3
       index=find(network_degreeFreq(:,2)<1.5);
       network_degreeFreq(index,2)=network_degreeFreq(index,2)+1;
    end
    [P,S]=polyfit(network_degreeFreq(:,1)',network_degreeFreq(:,2)',1);
    [y,DELTA]=polyval(P,network_degreeFreq(:,1)',S);
    figure(i);
    para=plot(network_degreeFreq(:,1)',y,'--ro',network_degreeFreq(:,1)',network_degreeFreq(:,2)','r*');
%     legend(strcat('paw low fit line and lamda is ',num2str(abs(P(1)))),'degree ditribution');
    set(para(1),'LineWidth',2);
    hold on;
    
    tempdata=network_statisticTable{GroupreactionDegree(i)};
    cluster_coefficient_groupreaction{i}=tempdata(:,2);
    network_degreeFreq=tabulate(tempdata(:,5));
    temp=[];
    index=1;
    for j=1:size(network_degreeFreq,1)
        if network_degreeFreq(j,1)~=0
            temp(index,:)=network_degreeFreq(j,:);
            index=index+1;
        end
    end
    network_degreeFreq=temp;
    network_degreeFreq=log(network_degreeFreq(:,[1,2]));
    degreecutoff=find(network_degreeFreq(:,1)<=3.4);
    network_degreeFreq=network_degreeFreq(degreecutoff,:);
    [P1,S1]=polyfit(network_degreeFreq(:,1)',network_degreeFreq(:,2)',1);
    [y,DELTA]=polyval(P1,network_degreeFreq(:,1)',S1);
    para1=plot(network_degreeFreq(:,1)',y,'--bd',network_degreeFreq(:,1)',network_degreeFreq(:,2)','b^');
    set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
    str1_legned=strcat({'\gamma:'},{num2str(abs(dataorgin(i,9)))});
    str2_legned=strcat({'\gamma:'},{num2str(abs(dataorgin(i,8)))});
    h_legned=legend([para(1),para1(1)],str1_legned{1},str2_legned{1});
    set(h_legned,'Fontsize',15,'FontWeight','bold');
    legend boxoff;
    set(para1(1),'LineWidth',2);
    set(para1(2),'LineWidth',2);
%     xlabel('The log scale of index of node degree','Fontsize',15,'FontWeight','bold');
%     ylabel('The frequency of node degree','Fontsize',15,'FontWeight','bold');
    species_name=filenameList(cofactorFreeDegree(i)+2).name;
    species_name=regexp(species_name,'\.','split');
%     title(strcat('The power law distribution for ',species_name{1}),'Fontsize',15,'FontWeight','bold');
    saveas(gcf,strcat('D:\huansheng\program\biggstat\fomulagraph\supplementary figure\20170425\',num2str(i),'.pdf'));
    saveas(gcf,strcat('D:\huansheng\program\biggstat\fomulagraph\supplementary figure\20170425\',num2str(i),'.fig'));
    saveas(gcf,strcat('D:\huansheng\program\biggstat\fomulagraph\supplementary figure\20170425\',num2str(i),'.emf'));
end

figure(19)
[P1,S1]=polyfit(network_degreeFreq(:,1)',network_degreeFreq(:,2)',1);
[y,DELTA]=polyval(P1,network_degreeFreq(:,1)',S1);
para1=plot(network_degreeFreq(:,1)',y,'--bd',network_degreeFreq(:,1)',network_degreeFreq(:,2)','b^');
hold on;
[P1,S1]=polyfit(network_degreeFreq(:,1)',2*network_degreeFreq(:,2)',1);
[y,DELTA]=polyval(P1,network_degreeFreq(:,1)',S1);
para1=plot(network_degreeFreq(:,1)',y,'--rd',network_degreeFreq(:,1)',2*network_degreeFreq(:,2)','r^');



[pvalue,hvalue]=ttest(dataorgin(:,8),dataorgin(:,9));


dataregress=dataorgin(:,[2,3,5,6,]);
% dataregress=dataorgin(:,[3,6,]);
dataregress=[dataregress';ones(1,size(dataregress,1))];
dataregress=dataregress';
dataregress(:,1)=dataregress(:,1)/1000;
dataregress(:,3)=dataregress(:,3)/1000;
datarestult=dataorgin(:,7);
datarestult=datarestult+0.25;

[b, bint,r,rint,stats]=regress(datarestult,dataregress(:,1:4));
[XL,YL,XS,YS,BATA,PCTVAT,MSE,STAT]=plsregress(datarestult,dataregress(:,1:4));

X1=rand(1,50)*0.3+0.6;
X2=rand(1,50)*0.3+0.1;
X3=ones(1,50);
X=[X1;X2;X3];
X=X';
% Y=X*b;
% figure(20);
% plot3(dataregress(:,1)',dataregress(:,2)',datarestult,'b*',X(:,1)',X(:,2),Y','r');
[B,I]=sort(dataorgin(:,8));
tempdis=dataorgin(:,9);

mean(tempdis)
mean(B)
powerlawcoefficient=[tempdis',B'];
powerlawproperty={};
index=1;

for i=1:length(tempdis)
    powerlawproperty{index}='CFF network';
    index=index+1;
end
for i=1:length(B)
    powerlawproperty{index}='PRF network';
    index=index+1;
end
figure(25)
para=boxplot(powerlawcoefficient,powerlawproperty);
set(gca,'Fontsize',15,'FontWeight','bold','linew',2);
set(gca,'xticklabel',{'1' '3'});
ylabel('Power law coefficient','Fontsize',15,'FontWeight','bold');
% xlabel('CFF network                     PRF network','Fontsize',15,'FontWeight','bold');
xlabel('CFF                            SRN','Fontsize',15,'FontWeight','bold');
set(para(:,:),'LineWidth',2);
set(findobj(gca,'Type','text'),'FontSize',15,'FontWeight','bold');



figure(21)
% axes('position',[0.68,0.64,0.2,0.25]);

x_data=dataregress(:,4);
x_data=x_data(I);
tempdis=tempdis(I);
% para=plot(1:length(B),B,'sr',1:length(B),tempdis,'sb',1:length(B),ones(1,length(B))*2,'-.k');
% para=plot(x_data,B,'sr',x_data,tempdis,'sb',min(x_data):(max(x_data)-min(x_data))/10:max(x_data),ones(1,11)*2,'-.k');
para=plot(x_data,B,'sr',x_data,tempdis,'sb');
Corr1=corr(x_data,B);
Corr2=corr(x_data,tempdis);
Corr1=corr(x_data,B);
Corr2=corr(x_data,tempdis);
Corr1=round(Corr1*10^3);
Corr1=Corr1/10^3;
Corr2=round(Corr2*10^3);
Corr2=Corr2/10^3;
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
set(para(1),'LineWidth',1.5);
set(para(2),'LineWidth',1.5);
h_legned=legend([para(1),para(2)],strcat('R=',num2str(Corr1)),strcat('R=',num2str(Corr2)));
set(h_legned,'Fontsize',15,'FontWeight','bold');
% legend boxoff;
xlabel('Proportion of metabolic genes','Fontsize',15,'FontWeight','bold');
ylabel('Power law coefficient','Fontsize',15,'FontWeight','bold');
figure(20)
% axes('position',[0.44,0.64,0.2,0.25]);
x_data=dataregress(:,2);
x_data=x_data(I);
% para=plot(x_data,B,'sr',x_data,tempdis,'sb',min(x_data):(max(x_data)-min(x_data))/10:max(x_data),ones(1,11)*2,'-.k');
para=plot(x_data,B,'sr',x_data,tempdis,'sb');
Corr1=corr(x_data,B);
Corr2=corr(x_data,tempdis);
Corr1=round(Corr1*10^3);
Corr1=Corr1/10^3;
Corr2=round(Corr2*10^3);
Corr2=Corr2/10^3;
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
set(gca,'xdir','reverse');
set(para(1),'LineWidth',1.5);
set(para(2),'LineWidth',1.5);
% h_legned=legend([para(1),para(2)],'PRF network','CFF network');
% h_legned=legend([para(1),para(2)],strcat('Corr:',num2str(Corr1)),strcat('Corr:',num2str(Corr2)));
h_legned=legend([para(1),para(2)],strcat('R=',num2str(Corr1)),strcat('R=',num2str(Corr2)));
set(h_legned,'Fontsize',15,'FontWeight','bold');

% legend boxoff;
xlabel('Ratio of #metabolites to #reacions','Fontsize',15,'FontWeight','bold');
ylabel('Power law coefficient','Fontsize',15,'FontWeight','bold');
%% clustering coefficient
filenameList_orginalstat = dir('D:\huansheng\program\biggstat\originalnetwork');
orginalnetwork_statisticTable={};
index=1;
for i=3:size(filenameList_orginalstat,1)
    orginalnetwork_statisticTable{index} = csvread(strcat('D:\huansheng\program\biggstat\originalnetwork\',filenameList_orginalstat(i).name),1,0);
    index=index+1;
end
coustering_coefficient_average=[];
coustering_coefficient_average_property={};
index=1;
for i=1:length(orginalnetwork_statisticTable)
    datatemp=orginalnetwork_statisticTable{i};
    coustering_coefficient_average(index)=mean(datatemp(:,2));
    coustering_coefficient_average_property{index}='Original';
    index=index+1;
end
for i=1:length(cluster_coefficient_cofactorfree)
    datatemp=cluster_coefficient_cofactorfree{i};
    coustering_coefficient_average(index)=mean(datatemp);
    coustering_coefficient_average_property{index}='CFF';
    index=index+1;
end
for i=1:length(cluster_coefficient_groupreaction)
    datatemp=cluster_coefficient_groupreaction{i};
    coustering_coefficient_average(index)=mean(datatemp);
    coustering_coefficient_average_property{index}='PRF';
    index=index+1;
end
figure(23)
para=boxplot(coustering_coefficient_average,coustering_coefficient_average_property);
set(gca,'Fontsize',15,'FontWeight','bold','linew',2);
set(gca,'xticklabel',{'1' '3'});
xlabel(' Original               CFF                  SRN    ','Fontsize',15,'FontWeight','bold');
ylabel('Clustering coefficient','Fontsize',15,'FontWeight','bold');
set(para(:,:),'LineWidth',2);
set(findobj(gca,'Type','text'),'FontSize',15,'FontWeight','bold');



random_network_stat_data=xlsread('RandomGraph_static.xlsx','RandomGraph_static','B2:E2334');
clustering_coefficient_ecoli=[];
clustering_coefficient_ecoli_proterty={};
index=0;
for i=3:length(filenameList_orginalstat)
    str=regexp(filenameList_orginalstat(i).name,'ecoli','match');
    if ~isempty(str)
        index=i-2;
        break;
    end
end
tempdata=orginalnetwork_statisticTable{index};
tempdata=tempdata(:,2);
tempdata=tempdata';
clustering_coefficient_ecoli=[clustering_coefficient_ecoli,tempdata];
index_property=1;
for i=1:length(tempdata)
    clustering_coefficient_ecoli_proterty{index_property}='Original';
    index_property=index_property+1;
end
index=0;
for i=3:length(filenameList)
    str=regexp(filenameList(i).name,'ecoli','match');
    if ~isempty(str)
        index=i-2;
        break;
    end
end
index=ceil(index/2);
tempdata=cluster_coefficient_cofactorfree{index};
tempdata=tempdata';
clustering_coefficient_ecoli=[clustering_coefficient_ecoli,tempdata];
for i=1:length(tempdata)
    clustering_coefficient_ecoli_proterty{index_property}='CFF';
    index_property=index_property+1;
end
tempdata=cluster_coefficient_groupreaction{index};
cluster_coefficient_groupreaction_data=tempdata;
tempdata=tempdata';
clustering_coefficient_ecoli=[clustering_coefficient_ecoli,tempdata];
for i=1:length(tempdata)
    clustering_coefficient_ecoli_proterty{index_property}='PRF';
    index_property=index_property+1;
end
tempdata=random_network_stat_data(:,2);
tempdata=tempdata';
clustering_coefficient_ecoli=[clustering_coefficient_ecoli,tempdata];
for i=1:length(tempdata)
    clustering_coefficient_ecoli_proterty{index_property}='BA';
    index_property=index_property+1;
end
figure(24);
% axes('position',[0.5,0.62,0.4,0.3]);
para=boxplot(clustering_coefficient_ecoli,clustering_coefficient_ecoli_proterty);
set(gca,'Fontsize',15,'FontWeight','bold','linew',2);
set(para,'LineWidth',2);
set(gca,'xticklabel',{'1' '3'});
xlabel('      Original         CFF           SRN             BA         ','Fontsize',15,'FontWeight','bold');
ylabel('Clustering coefficient','Fontsize',15,'FontWeight','bold');
set(findobj(gca,'Type','text'),'FontSize',15,'FontWeight','bold');

figure(27)
% [P1,S1]=polyfit(dataorgin(:,3),dataorgin(:,10),1);
% [y1,DELTA]=polyval(P1,dataorgin(:,3),S1);
% [P2,S2]=polyfit(dataorgin(:,3),dataorgin(:,11),1);
% [y2,DELTA]=polyval(P2,dataorgin(:,3),S2);
% para=plot(dataorgin(:,3),dataorgin(:,10),'sr',dataorgin(:,3),dataorgin(:,11),'sb',dataorgin(:,3),y1,'--r',dataorgin(:,3),y2,'--b');
para=plot(dataorgin(:,3),dataorgin(:,10),'sr',dataorgin(:,3),dataorgin(:,11),'sb');
Corr1=corr(dataorgin(:,3),dataorgin(:,10));
Corr2=corr(dataorgin(:,3),dataorgin(:,11));
Corr1=round(Corr1*10^3);
Corr1=Corr1/10^3;
Corr2=round(Corr2*10^3);
Corr2=Corr2/10^3;
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
set(gca,'xdir','reverse')
set(para,'LineWidth',1.5);
% h_legned=legend(para,'PRF network');
% h_legned=legend(para,strcat('PRF network,Corr= ',num2str(Corr1)));
% h_legned=legend([para(1),para(2)],strcat('Corr:',num2str(Corr1)),strcat('Corr:',num2str(Corr2)));
h_legned=legend([para(1),para(2)],strcat('R=',num2str(Corr1)),strcat('R=',num2str(Corr2)));

set(h_legned,'Fontsize',15,'FontWeight','bold');
xlabel('Ratio of #metabolites to #reacions','Fontsize',15,'FontWeight','bold');
ylabel('Clustering coefficient','Fontsize',15,'FontWeight','bold');

figure(28)
para=plot(dataorgin(:,6),dataorgin(:,10),'sr',dataorgin(:,6),dataorgin(:,11),'sb');
Corr1=corr(dataorgin(:,6),dataorgin(:,10));
Corr2=corr(dataorgin(:,6),dataorgin(:,11));
Corr1=round(Corr1*10^3);
Corr1=Corr1/10^3;
Corr2=round(Corr2*10^3);
Corr2=Corr2/10^3;
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
set(para,'LineWidth',1.5);
% h_legned=legend(para,'PRF network');
% h_legned=legend(para,strcat('PRF network,Corr= ',num2str(Corr1)));
% h_legned=legend([para(1),para(2)],strcat('Corr:',num2str(Corr1)),strcat('Corr:',num2str(Corr2)));
h_legned=legend([para(1),para(2)],strcat('R=',num2str(Corr1)),strcat('R=',num2str(Corr2)));
set(h_legned,'Fontsize',15,'FontWeight','bold');
xlabel('Proportion of metabolic genes','Fontsize',15,'FontWeight','bold');
ylabel('Clustering coefficient','Fontsize',15,'FontWeight','bold');


figure(22)
CCSimulate=load('clustering coefficient simulate.txt');
para=plot(CCSimulate(:,1),CCSimulate(:,3),'dr',CCSimulate(:,1),CCSimulate(:,3),'r',CCSimulate(:,1),CCSimulate(:,4),'*b',CCSimulate(:,1),CCSimulate(:,4),'b');
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
set(para(1),'LineWidth',1.5);
set(para(2),'LineWidth',1.5);
set(para(3),'LineWidth',1.5);
set(para(4),'LineWidth',1.5);
h_legend=legend([para(2),para(4)],'Average degree of neighboring node','Interconnected edges of neighboring node');
set(h_legned,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold');


% data_x=1:200;
data_x=1:50;
data_y(1:50)=data_x(1:50)*0.01;
% data_y(51:200)=0.5-(data_x(51:200)-51)*0.0033;
axes('position',[0.7,0.16,0.2,0.29]);
para=plot(log(data_x),data_y,'*r');
set(gca,'Fontsize',13,'FontWeight','bold','linew',1.5);
set(para,'LineWidth',1.5);

degree_groupreaction_data=network_statisticTable{8};
% degree_groupreaction_data=network_statisticTable{6};
cluster_coefficient_groupreaction_data=degree_groupreaction_data(:,2);
degree_groupreaction_data=degree_groupreaction_data(:,5);
max_index=max(degree_groupreaction_data);
clusteringcoef_degree=[];
for i=0:max_index
    find_label=find(degree_groupreaction_data==i);
    if isempty(find_label)
        continue;
    end
    temp=sum(cluster_coefficient_groupreaction_data(find_label))/length(find_label);
    clusteringcoef_degree=[clusteringcoef_degree;i,temp];
end
figure(30)
para=plot(log(clusteringcoef_degree(:,1)),clusteringcoef_degree(:,2),'r*');
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
xlabel('Log-transformed node-degree','FontSize',15,'FontWeight','bold');
ylabel('Average clustering coefficient','FontSize',15,'FontWeight','bold');
% title('E coli');
set(para,'LineWidth',1.5);
corr(log(clusteringcoef_degree(2:end,1)),clusteringcoef_degree(2:end,2))
%% the neighbor degree distributin for real network
data=load('reaction_group_graph_STM_v1_0.txt');
% data=load('reaction_group_graph_iAF692.txt');
max_index=max(data);
max_index=max(max_index)
min_index=min(data);
min_index=min(min_index);

outputDegree=[];
neighbor={};
for i=min_index:max_index
    degree_index=find(data(:,1)==i);
    if length(degree_index)>0
        outputDegree(i)=length(degree_index);
        neighbor{i}=data(degree_index,2);
    else
        outputDegree(i)=0;
        neighbor{i}=0;
    end
end

max_index=max(outputDegree);
min_index=min(outputDegree);
neighbor_degree_static=[];
for i=1:max_index
    temp=find(outputDegree==i);
    if length(temp)>0
        neighbor_temp=[];
        for j=1:length(temp)
            neighbor_temp=[neighbor_temp,neighbor{temp(j)}'];
        end
        neighbor_temp=unique(neighbor_temp);
        neighbor_temp_vaerage=mean(outputDegree(neighbor_temp));
        neighbor_degree_static(i,1)=i;
        neighbor_degree_static(i,2)=neighbor_temp_vaerage;
    else
%         neighbor_degree_static(i,1)=i;
%         neighbor_degree_static(i,2)=0;
    end
end
% ref_degree=sqrt(2*(clusteringcoef_degree(:,2))).*clusteringcoef_degree(:,1);
ref_degree=(2*(clusteringcoef_degree(:,2))).*clusteringcoef_degree(:,1);
ref_degree=sqrt(2.2*(clusteringcoef_degree(:,1)-1).*(clusteringcoef_degree(:,2).*clusteringcoef_degree(:,2).*(clusteringcoef_degree(:,1)-1)/2+(1-clusteringcoef_degree(:,2)/2)));
figure(100)
para1=plot(neighbor_degree_static(:,1),neighbor_degree_static(:,2),'rd');
hold on;
para2=plot(clusteringcoef_degree(:,1),ref_degree,'b*');
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
xlabel('The degree of focused node');
ylabel('The average degree of neighboring node');
% title('E coli');
set(para1,'LineWidth',2);
set(para2,'LineWidth',2);
legend([para1,para2],'Real degree of neighboring node','Predicted degree of neighboring node');
% legend(para2,'predict average neighbor degree of E coli');