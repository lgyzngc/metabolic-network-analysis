% cofactor={'adp[c]','atp[c]','amp[c]','nadp[c]','nadph[c]','nad[c]','nadh[c]','h2o[c]','h[c]','h[e]','coa[c]','ACP[c]','pi[c]','ppi[c]','cmp[c]','q8h2[c]','q8[c]','mqn8[c]','mql8[c]','gtp[c]','gdp[c]','ctp[c]','pi[e]','udp[c]','fadh2[c]','gmp[c]','ump[c]'};
% cofactor={'adp[c]','atp[c]','amp[c]','nadp[c]','nadph[c]','nad[c]','nadh[c]','h2o[c]','h2o[p]','h[c]','h[e]','h[p]','pi[c]','pi[p]','ppi[c]','cmp[c]','gtp[c]','gdp[c]','ctp[c]','pi[e]','udp[c]','fadh2[c]','gmp[c]','ump[c]','fe2[c]'};
cofactor={'adp[c]','atp[c]','amp[c]','nadp[c]','nadph[c]','nad[c]','nadh[c]','h2o[c]','h2o[p]','h[c]','h[e]','h[p]','coa[c]','ACP[c]','pi[c]','pi[p]','ppi[c]','cmp[c]','q8h2[c]','q8[c]','mqn8[c]','mql8[c]','gtp[c]','gdp[c]','ctp[c]','pi[e]','udp[c]','fadh2[c]','gmp[c]','ump[c]','fe2[c]','co2[c]','co2[e]','co2[p]',...
    'o2[c]','adp[g]','atp[g]','amp[g]','nadp[g]','nadph[g]','nad[g]','nadh[g]','h2o[g]','h[g]','pi[g]','gtp[g]','gdp[g]','ctp[g]','udp[g]','fadh2[g]','gmp[g]','ump[g]','co2[g]','o2[g]'...
    ,'adp[m]','atp[m]','amp[m]','nadp[m]','nadph[m]','nad[m]','nadh[m]','h2o[m]','h[m]','pi[m]','gtp[m]','gdp[m]','ctp[m]','udp[m]','fadh2[m]','gmp[m]','ump[m]','co2[m]','o2[m]'...
    ,'adp[n]','atp[n]','amp[n]','nadp[n]','nadph[n]','nad[n]','nadh[n]','h2o[n]','h[n]','pi[n]','gtp[n]','gdp[n]','ctp[n]','udp[n]','fadh2[n]','gmp[n]','ump[n]','co2[n]','o2[n]'...
    ,'adp[v]','atp[v]','amp[v]','nadp[v]','nadph[v]','nad[v]','nadh[v]','h2o[v]','h[v]','pi[v]','gtp[v]','gdp[v]','ctp[v]','udp[v]','fadh2[v]','gmp[v]','ump[v]','co2[v]','o2[v]'...
    ,'adp[x]','atp[x]','amp[x]','nadp[x]','nadph[x]','nad[x]','nadh[x]','h2o[x]','h[x]','pi[x]','gtp[x]','gdp[x]','ctp[x]','udp[x]','fadh2[x]','gmp[x]','ump[x]','co2[x]','o2[x]'...
    };
% cofactor={'q8h2[c]'};
import java.util.Hashtable;
% Ecoli_Model = readCbModel('Bs_iYO844_flux1.xml');
% Ecoli_Model = readCbModel('iIT341_min_medI.xml');
% Ecoli_Model = readCbModel('Mb_iAF692.xml');
% Ecoli_Model = readCbModel('Ec_iAF1260_flux1.xml');
% Ecoli_Model = readCbModel('iNJ661_middlebrook.xml');


% Ecoli_Model = readCbModel('C:\Users\lyg\Dropbox\Network_opinion\xml_models\sbml3\iJO1366_2.xml');
% Ecoli_Model = readCbModel('C:\Users\lyg\Dropbox\Network_opinion\xml_models\sbml\iMM904.xml');

% Ecoli_Model=model;
cofactor_table = java.util.Hashtable;
for i=1:length(cofactor)
    cofactor_table.put(cofactor{i},i);
end
% load iHN637.mat;
% load iMM904.mat;
load D:\huansheng\program\bigg\iAF987.mat;
 

Ecoli_Model=iAF987;

% Ecoli_Model = readCbModel('Ec_iAF1260_flux1.xml');

for i=1:length(Ecoli_Model.mets)
    a=regexp(Ecoli_Model.mets{i},'_[cempgnvx]$','match');
    if ~isempty(a)
        Ecoli_Model.mets{i}=regexprep(Ecoli_Model.mets{i},'_[cempgnvx]$','[$0]');
        Ecoli_Model.mets{i}=strrep(Ecoli_Model.mets{i},'[_','[');
        disp(Ecoli_Model.mets{i})
    end
end
for i=1:length(Ecoli_Model.grRules)
    a=regexp(Ecoli_Model.grRules{i},'_','match');
    if ~isempty(a)
        Ecoli_Model.grRules{i}=regexprep(Ecoli_Model.grRules{i},'_','csbl');
    end
end


network_table={};
input_network_table={};
output_network_table={};
Reaction_degree=[];
Reaction_input_degree=[];
Reaction_output_degree=[];
reac_index=[];
index=1;
in_index=1;
out_index=1;
for i=1:size(Ecoli_Model.S,2)
    index_set=find(Ecoli_Model.S(:,i)~=0);
    meta_index_set=[];
    meta_input_index=[];
    meta_output_index=[];
    
    meta_input_set=[];
    meta_output_set=[];
    reac_index=[];
    for j=1:length(index_set)
        if ~cofactor_table.containsKey(Ecoli_Model.mets{index_set(j)})
            if Ecoli_Model.rev(i)==1
                
%                 meta_index_set=setdiff(meta_index_set,i);
                meta_index_set=find(Ecoli_Model.S(index_set(j),:)<0);
                temp_reverse=find(Ecoli_Model.S(index_set(j),:)>0);
                temp_reverse=temp_reverse(find(Ecoli_Model.rev(temp_reverse)==1));
                if ~isempty(temp_reverse)
                    meta_index_set=[meta_index_set,temp_reverse];
                end
                meta_output_index=meta_index_set;
                meta_output_set=[meta_output_set,meta_output_index];
                
                meta_index_set=find(Ecoli_Model.S(index_set(j),:)>0);
                temp_reverse=find(Ecoli_Model.S(index_set(j),:)<0);
                temp_reverse=temp_reverse(find(Ecoli_Model.rev(temp_reverse)==1));
                if ~isempty(temp_reverse)
                    meta_index_set=[meta_index_set,temp_reverse];
                end
                meta_input_index=meta_index_set;
                meta_input_set=[meta_input_set,meta_input_index];
            elseif Ecoli_Model.S(index_set(j),i)>0
                meta_index_set=find(Ecoli_Model.S(index_set(j),:)<0);
                temp_reverse=find(Ecoli_Model.S(index_set(j),:)>0);
                temp_reverse=temp_reverse(find(Ecoli_Model.rev(temp_reverse)==1));
                if ~isempty(temp_reverse)
                    meta_index_set=[meta_index_set,temp_reverse];
                end
                meta_output_index=meta_index_set;
                meta_output_set=[meta_output_set,meta_output_index];
            elseif Ecoli_Model.S(index_set(j),i)<0
                meta_index_set=find(Ecoli_Model.S(index_set(j),:)>0);
                temp_reverse=find(Ecoli_Model.S(index_set(j),:)<0);
                temp_reverse=temp_reverse(find(Ecoli_Model.rev(temp_reverse)==1));
                if ~isempty(temp_reverse)
                    meta_index_set=[meta_index_set,temp_reverse];
                end
                meta_input_index=meta_index_set;
                meta_input_set=[meta_input_set,meta_input_index];
            end
            reac_index=[reac_index,meta_index_set];
%             meta_input_set=[meta_input_set,meta_input_index];
%             meta_output_set=[meta_output_set,meta_output_index];
        else
            cofactor_table.put(Ecoli_Model.mets{index_set(j)},index_set(j));
        end
    end
    
    
    
    reac_index=unique(reac_index);
    repeat_i=find(reac_index==i);
    
    meta_input_set=unique(meta_input_set);
    repeat_input_i=find(meta_input_set==i);
    if ~isempty(repeat_input_i)
        Reaction_input_degree(i)=length(meta_input_set)-1;    
        meta_input_set(repeat_input_i)=[];
    else
        Reaction_input_degree(i)=length(meta_input_set);
    end
    
    meta_output_set=unique(meta_output_set);
    repeat_output_i=find(meta_output_set==i);
    if ~isempty(repeat_output_i)
        Reaction_output_degree(i)=length(meta_output_set)-1;    
        meta_output_set(repeat_output_i)=[];
    else
        Reaction_output_degree(i)=length(meta_output_set);
    end
    
    
    if ~isempty(repeat_i)
        Reaction_degree(i)=length(reac_index)-1;    
        reac_index(repeat_i)=[];
    else
        Reaction_degree(i)=length(reac_index);
    end
%     for j=1:length(reac_index)
%         if cofactor_table.containsKey(Ecoli_Model.mets{reac_index(j)})
%             reac_index(j)=[];
%         end
%     end
    for j=1:length(reac_index)
        network_table{index,1}=Ecoli_Model.rxns{i};
        network_table{index,2}=Ecoli_Model.rxns{reac_index(j)};
        network_table{index,3}='1';
        index=index+1;
    end
    

    for j=1:length(meta_input_set)
%         input_network_table{in_index,1}=Ecoli_Model.rxns{meta_input_set(j)};
%         input_network_table{in_index,2}=Ecoli_Model.rxns{i};
        input_network_table{in_index,1}=meta_input_set(j);
        input_network_table{in_index,2}=i;
        input_network_table{in_index,3}='1';
        in_index=in_index+1;
    end
    
    for j=1:length(meta_output_set)
        
%         output_network_table{out_index,1}=Ecoli_Model.rxns{i};
%         output_network_table{out_index,2}=Ecoli_Model.rxns{meta_output_set(j)};
        output_network_table{out_index,1}=i;
        output_network_table{out_index,2}=meta_output_set(j);
        output_network_table{out_index,3}='1';
        out_index=out_index+1;
    end
end
graph_output=fopen('network_without_cofactor.txt','w');
for i=1:size(input_network_table,1)
    fprintf(graph_output,'%d\t%d\n',input_network_table{i,1},input_network_table{i,2});
end
for i=1:size(output_network_table,1)
    fprintf(graph_output,'%d\t%d\n',output_network_table{i,1},output_network_table{i,2});
end
fclose(graph_output);


figure(1)
bar(sort(Reaction_degree));

max_index=max(Reaction_degree);
reaction_degree_dis=[];
for i=1:max_index
    reaction_degree_dis(i)=length(find(Reaction_degree == i));
end
figure(2)
plot(1:max_index,reaction_degree_dis,'*');

file_ID=fopen('reaction_network.txt','w');
for i=1:size(network_table,1)
    fprintf(file_ID,'%d\t%d\t%s\n',network_table{i,1},network_table{i,2},network_table{i,3});
end
fclose(file_ID);

file_ID=fopen('input_reaction_network_iSB619.txt','w');
for i=1:size(input_network_table,1)
    fprintf(file_ID,'%d\t%d\t%s\n',input_network_table{i,1},input_network_table{i,2},input_network_table{i,3});
end
fclose(file_ID);
file_ID=fopen('input_reaction_network_network.txt','w');
for i=1:size(input_network_table,1)
    fprintf(file_ID,'%s\t%s\t%s\n',Ecoli_Model.rxns{input_network_table{i,1}},Ecoli_Model.rxns{input_network_table{i,2}},input_network_table{i,3});
end
fclose(file_ID);

file_ID=fopen('output_reaction_network.txt','w');
for i=1:size(output_network_table,1)
    fprintf(file_ID,'%d\t%d\t%s\n',output_network_table{i,1},output_network_table{i,2},output_network_table{i,3});
end
fclose(file_ID);
%% enzyme dis (and or single) and degree dis, and for muti-function raction, or for evolutionary in different conditions, single gene for enssentile.
% single with 1, and with 2, or with 3, and_or with 4
enzyme_single_hash = java.util.Hashtable;
enzyme_and_hash = java.util.Hashtable;
enzyme_or_hash = java.util.Hashtable;
enzyme_and_or_hash = java.util.Hashtable;
enzyme_index_hash = java.util.Hashtable;
hash_key_index={};


enzyme_cluster_set=[];
enzyme_cluster_stic=[];
enzyme_reaction_match=[];
enzyme_index=1;
for i=1:length(Ecoli_Model.grRules)
    S = regexp(Ecoli_Model.grRules{i}, 'or', 'split');
    if isempty(S{1})
        continue;
    end
    
%             pattern_gene_temp=regexp(Ecoli_Model.grRules{i},'\w(\d{4})','match');%parrern for Ecoli gene ID;
            pattern_gene_temp=regexp(Ecoli_Model.grRules{i},'(\w+)(\d+)','match');%parrern for other gene ID;
%            pattern_gene_temp=regexp(Ecoli_Model.grRules{i},'(\w)+(\d)+','match');
%              pattern_gene_temp=regexp(Ecoli_Model.grRules{i},'(\<Rv[\d\w]+\>)|((\d)+\.(\d)+)','match');%for human 
            if isempty(pattern_gene_temp)
                continue;
            end
            
            pattern_gene_temp=sort(pattern_gene_temp);
            pattern_gene_set='';
            for pg=1:length(pattern_gene_temp)
                pattern_gene_set=strcat(pattern_gene_set,'_',pattern_gene_temp{pg});
            end
            hash_key_index{i}=pattern_gene_set;
            if isEmpty(enzyme_index_hash)
                enzyme_index_hash.put(pattern_gene_set,enzyme_index);
                enzyme_reaction_match(i)=enzyme_index;
                enzyme_index=enzyme_index+1;
                
            else
                if enzyme_index_hash.containsKey(pattern_gene_set)
                    enzyme_reaction_match(i)=enzyme_index_hash.get(pattern_gene_set);
                else
                    enzyme_index_hash.put(pattern_gene_set,enzyme_index);
                    enzyme_reaction_match(i)=enzyme_index;
                    enzyme_index=enzyme_index+1;
                end
            end
    
    if length(S) == 1 
        S_sub=regexp(S{1}, 'and', 'split');
        if length(S_sub) == 1
            enzyme_cluster_set(i)=1;
            enzyme_cluster_stic(i)=1;
            
            if isEmpty(enzyme_single_hash)
                enzyme_single_hash.put(pattern_gene_set,1);
            else
                if enzyme_single_hash.containsKey(pattern_gene_set)
                    enzyme_single_hash.put(pattern_gene_set,1+enzyme_single_hash.get(pattern_gene_set));
                else
                    enzyme_single_hash.put(pattern_gene_set,1);
                end
            end
            
        else
            enzyme_cluster_set(i)=2;
            enzyme_cluster_stic(i)=length(S_sub);
            
            if isEmpty(enzyme_and_hash)
                enzyme_and_hash.put(pattern_gene_set,1);
            else
                if enzyme_and_hash.containsKey(pattern_gene_set)
                    enzyme_and_hash.put(pattern_gene_set,1+enzyme_and_hash.get(pattern_gene_set));
                else
                    enzyme_and_hash.put(pattern_gene_set,1);
                end
            end
            
        end
    else
        flag=0;
        for j=1:length(S)
            S_sub=regexp(S{j}, 'and', 'split');
            if length(S_sub)>1
                enzyme_cluster_set(i)=4;
                enzyme_cluster_stic(i)=30+length(S)+5*length(S_sub);
                flag=1;
                
            if isEmpty(enzyme_and_or_hash)
                enzyme_and_or_hash.put(pattern_gene_set,1);
            else
                if enzyme_and_or_hash.containsKey(pattern_gene_set)
                    enzyme_and_or_hash.put(pattern_gene_set,1+enzyme_and_or_hash.get(pattern_gene_set));
                else
                    enzyme_and_or_hash.put(pattern_gene_set,1);
                end
            end
                
                break;
            end
        end
        if flag==0
            enzyme_cluster_set(i)=3;
            enzyme_cluster_stic(i)=20+ length(S);
            
            if isEmpty(enzyme_or_hash)
                enzyme_or_hash.put(pattern_gene_set,1);
            else
                if enzyme_or_hash.containsKey(pattern_gene_set)
                    enzyme_or_hash.put(pattern_gene_set,1+enzyme_or_hash.get(pattern_gene_set));
                else
                    enzyme_or_hash.put(pattern_gene_set,1);
                end
            end
            
        end
    end
end
 if i<=length(Reaction_degree)
     for i=i:length(Reaction_degree)
         enzyme_cluster_set(i)=0;
         enzyme_cluster_stic(i)=0;
     end
 end


enzyme_single_set=enzyme_single_hash.values.toArray;
enzyme_and_set=enzyme_and_hash.values.toArray;
enzyme_or_set=enzyme_or_hash.values.toArray;
enzyme_and_or_set=enzyme_and_or_hash.values.toArray;

figure(9)
plot(Reaction_degree,enzyme_cluster_stic,'*');

figure(17)
whole_genes_num=[];
for i=1:length(hash_key_index)
    temp_name=hash_key_index{i};
    if isempty(temp_name)
        whole_genes_num(i)=0;
        continue;
    end
    if enzyme_single_hash.containsKey(temp_name)
        whole_genes_num(i)=enzyme_single_hash.get(temp_name);
    elseif enzyme_and_hash.containsKey(temp_name)
        whole_genes_num(i)=enzyme_and_hash.get(temp_name);    
    elseif enzyme_or_hash.containsKey(temp_name)
        whole_genes_num(i)=enzyme_or_hash.get(temp_name);
    elseif enzyme_and_or_hash.containsKey(temp_name)
        whole_genes_num(i)=enzyme_and_or_hash.get(temp_name);
    else
        whole_genes_num(i)=0;
    end
end
if length(whole_genes_num)<length(Reaction_degree)
    for i=length(whole_genes_num)+1:length(Reaction_degree)
        whole_genes_num(i)=0;
    end
end

plot(whole_genes_num,Reaction_degree,'*');



figure(24)
plot(1:length(whole_genes_num),whole_genes_num,'*');
figure(26)
plot(1:length(Reaction_degree),Reaction_degree,'*');

irres=find(Ecoli_Model.rev==0);
% irres=find(Ecoli_Model.rev>0);
B=enzyme_cluster_stic(irres);
[B,I]=sort(B);

C=whole_genes_num(irres);
C=C(I);

D=Reaction_degree(irres);
D=D(I);
figure(27);
plot(1:length(I),B,'*',1:length(I),C,'r^',1:length(I),D,'*g');

E=C(1474:length(C));
[E,IE]=sort(E);
F=D(1474:length(D));
F=F(IE);
figure(28)
boxplot(F,E);

%% enzyme gene num distribution and homology raction num
gene_vs_homo_reac=[];
index=1;
for i=1:length(B)
    if(C(i)>1)
        gene_vs_homo_reac(1,index)=C(i);
        temp=B(i);
        if temp==1
            gene_vs_homo_reac(2,index)=1;
        elseif temp>1 && temp<20
            gene_vs_homo_reac(2,index)=3;
        elseif temp<40
            gene_vs_homo_reac(2,index)=2;
        else
            gene_vs_homo_reac(2,index)=4;
        end
        index=index+1;
    end
end
for i=1:4
    temp_id=find(gene_vs_homo_reac(2,:)==i);
    if isempty(temp_id)
        gene_vs_homo_reac=[gene_vs_homo_reac';1,i];
        gene_vs_homo_reac=gene_vs_homo_reac';
    end
end
figure(29)
bh=boxplot(gene_vs_homo_reac(1,:)',gene_vs_homo_reac(2,:)','datalim',[0,55],'extrememode','compress' ,'notch','on','label',{'GE','MHGE','MSGE','MHSGE'});
set(gca,'Fontsize',13,'FontWeight','bold','linew',2);
ylabel('The distribution of homolog reaction','Fontsize',15,'FontWeight','bold');
set(bh(:,:),'LineWidth',2,'Color','k');
set(findobj(gca,'Type','text'),'FontSize',15,'FontWeight','bold');
% title('The enzyme complexity vs the distribution of homolog reaction');

gene_vs_homo_single_index=find(gene_vs_homo_reac(2,:)<3);
gene_vs_homo_muti_index=find(gene_vs_homo_reac(2,:)>2);
gene_vs_homo_reac_set=gene_vs_homo_reac(1,:);
[h_value,p_value]=ttest2(gene_vs_homo_reac_set(gene_vs_homo_single_index),gene_vs_homo_reac_set(gene_vs_homo_muti_index));
mean(gene_vs_homo_reac_set(gene_vs_homo_single_index))
mean(gene_vs_homo_reac_set(gene_vs_homo_muti_index))




%% combine reactions set which are more than 7 catazying by same enzyme to a common group reaction.
sort_replace_set=enzyme_reaction_match;
group_raction_set={};
index=1;
for i=1:length(sort_replace_set)
    if sort_replace_set(i)<=0
        continue;
    end
    temp=find(sort_replace_set==sort_replace_set(i));
    if length(temp)>5 && length(temp)~=67
        temp_genename=Ecoli_Model.grRules{i};
        sort_replace_set(temp)=-1;
        ractant_set=[];
        product_set=[];
        for j=1:length(temp)
            ractant_index=find(Ecoli_Model.S(:,temp(j))<0);
            product_index=find(Ecoli_Model.S(:,temp(j))>0);
%             if Ecoli_Model.rev(temp(j))==1 && j==1
            if Ecoli_Model.rev(temp(j))==1
                ractant_index=[ractant_index;product_index];
                product_index=[product_index;ractant_index];
            end
            ractant_set=[ractant_set;ractant_index];
            product_set=[product_set;product_index];
        end
        ractant_set=unique(ractant_set);
        product_set=unique(product_set);
        ractant_without_cofactor=[];
        for j=1:length(ractant_set)
            if ~cofactor_table.containsValue(ractant_set(j))
                ractant_without_cofactor=[ractant_without_cofactor,ractant_set(j)];
            end
        end
        product_without_cofactor=[];
        for j=1:length(product_set)
            if ~cofactor_table.containsValue(product_set(j))
                product_without_cofactor=[product_without_cofactor,product_set(j)];
            end
        end
        group_raction_set{index,1}=ractant_without_cofactor;
        group_raction_set{index,2}=product_without_cofactor;
        group_raction_set{index,3}=temp;
        group_raction_set{index,4}=length(temp);
        group_raction_set{index,5}={Ecoli_Model.metNames{ractant_without_cofactor}};
        group_raction_set{index,6}={Ecoli_Model.metNames{product_without_cofactor}};
        group_raction_set{index,7}=temp_genename;
        group_raction_set{index,8}=intersect(ractant_without_cofactor,product_without_cofactor);
        index=index+1;
    end
end
group_reaction_set_index=[];
for i=1:size(group_raction_set,1)
    group_reaction_set_index=[group_reaction_set_index,group_raction_set{i,3}];
end
%% group reaction output
group_reaction_file=fopen('group_reactions_ref_iMM904.txt','w');
for i=1:size(group_raction_set,1)
    substrates_index_set_groupreaction=group_raction_set{i,1};
    products_index_set_groupreaction=group_raction_set{i,2};
    groupreaction_index=group_raction_set{i,3};
    is_reversible=length(group_raction_set{i,8});
    fprintf(group_reaction_file,'groupreaction:\t');
    for j=1:length(substrates_index_set_groupreaction)
        if j==length(substrates_index_set_groupreaction)
            fprintf(group_reaction_file,'%s ',Ecoli_Model.metNames{substrates_index_set_groupreaction(j)});
        else
            fprintf(group_reaction_file,'%s + ',Ecoli_Model.metNames{substrates_index_set_groupreaction(j)});
        end
    end
    if is_reversible==1
        fprintf(group_reaction_file,'----> ');
    else
        fprintf(group_reaction_file,'<----> ');
    end
    for j=1:length(products_index_set_groupreaction)
        if j==length(products_index_set_groupreaction)
            fprintf(group_reaction_file,'%s ',Ecoli_Model.metNames{products_index_set_groupreaction(j)});
        else
            fprintf(group_reaction_file,'%s + ',Ecoli_Model.metNames{products_index_set_groupreaction(j)});            
        end
        fprintf(group_reaction_file,'%s + ',Ecoli_Model.metNames{products_index_set_groupreaction(j)});
    end
    fprintf(group_reaction_file,'\nsub reactions set:\t%s\n',group_raction_set{i,7});
    for j=1:length(groupreaction_index)
        substrate_indexset=find(Ecoli_Model.S(:,groupreaction_index(j))<0);
        product_indexset=find(Ecoli_Model.S(:,groupreaction_index(j))>0);
        if length(substrate_indexset)<1
            continue;
        end
        for k=1:length(substrate_indexset)
            if k==length(substrate_indexset)
                fprintf(group_reaction_file,'%s ',Ecoli_Model.metNames{substrate_indexset(k)});
            else
                fprintf(group_reaction_file,'%s + ',Ecoli_Model.metNames{substrate_indexset(k)});
            end
        end
        if Ecoli_Model.rev(substrate_indexset(k))==1
            fprintf(group_reaction_file,'----> ');
        else
            fprintf(group_reaction_file,'<----> ');
        end
        if length(product_indexset)<1
            continue;
        end
        for k=1:length(product_indexset)
            if k==length(product_indexset)
                fprintf(group_reaction_file,'%s ',Ecoli_Model.metNames{product_indexset(k)});
            else
                fprintf(group_reaction_file,'%s + ',Ecoli_Model.metNames{product_indexset(k)});
            end
        end
        fprintf(group_reaction_file,'\n');
    end
    fprintf(group_reaction_file,'\n');
    fprintf(group_reaction_file,'\n');
end
fclose(group_reaction_file);

%% reconstruct new network
recons_s=Ecoli_Model.S;
for i=1:length(cofactor)
    recons_s(cofactor_table.get(cofactor{i}),:)=0;
end
reaction_preperity=[];
for i=1:size(recons_s,2)
    reaction_preperity(i,1)=i;
    reaction_preperity(i,2)=0;
    meta_temp=find(recons_s(:,i)~=0);
    reaction_preperity(i,3)=length(meta_temp);
    reaction_preperity(i,4)=Ecoli_Model.rev(i);
    reaction_preperity(i,5)=0;
end
for i=1:size(group_raction_set,1)
    reaction_preperity(group_raction_set{i,3},5)=i;
end
file_preperty=fopen('reaction_preperty_iMM904.txt','w');
for i=1:size(reaction_preperity,1)
%     fprintf(file_preperty,'%d = %d\n',i,reaction_preperity(i,1));
    fprintf(file_preperty,'%d\t%d\t%d\t%d\t%d\n',reaction_preperity(i,1),reaction_preperity(i,2),reaction_preperity(i,3),reaction_preperity(i,4),reaction_preperity(i,5));
end
fclose(file_preperty);
meta_len=size(recons_s,1);
delete_Reactionindex=[];
for i=1:size(group_raction_set,1)
    temp=(group_raction_set(i,3));
    recons_s(:,temp{1})=0;
    delete_Reactionindex=[delete_Reactionindex,temp{1}];
    index_substrate=group_raction_set(i,1);
    index_product=group_raction_set(i,2);
    
    if isempty(intersect(index_substrate{1},index_product{1}))
        catabolism_array=zeros(meta_len,1);
        catabolism_array(index_substrate{1})=-1;
        catabolism_array(index_product{1})=1;
        group_raction_set{i,8}=size(recons_s,2)+1;
    else
        catabolism_array=zeros(meta_len,2);
  %      intersection=intersect(group_raction_set(i,1),group_raction_set(i,2));
        catabolism_array(index_product{1},1)=1;
        catabolism_array(index_substrate{1},1)=-1;
        
        
        catabolism_array(index_substrate{1},2)=1;
        catabolism_array(index_product{1},2)=-1;
        group_raction_set{i,8}=[size(recons_s,2)+1,size(recons_s,2)+2];
    end
    recons_s=[recons_s,catabolism_array];
end

subnetwork_properity=[];
for i=1:size(recons_s,2)
    subnetwork_properity(i)=0;
end

reconstrut_graph=[];
graph_index=1;
reac_len=size(recons_s,2);
for i=1:reac_len
    temp_p=find(recons_s(:,i)>0);
    temp_r=find(recons_s(:,i)<0);

    if i>length(Ecoli_Model.rev) || Ecoli_Model.rev(i)==0
        candidate_connect=[];
        for j=1:length(temp_p)
%         for j=1:length(temp_r)
%             candidate_connect=[candidate_connect,find(recons_s(temp_r(j),:)>0)];
%             temp_reverse=find(recons_s(temp_r(j),:)<0);
            
            candidate_connect=[candidate_connect,find(recons_s(temp_p(j),:)<0)];
            temp_reverse=find(recons_s(temp_p(j),:)>0);
             
            temp_reverse=temp_reverse(find(temp_reverse<=length(Ecoli_Model.rxns)));
            temp_reverse=temp_reverse(find(Ecoli_Model.rev(temp_reverse)==1));
            if  ~isempty(temp_reverse)
%                 candidate_connect=[candidate_connect,temp_reverse];
            end
        end
        
        candidate_connect=unique(candidate_connect);
        candidate_connect=setdiff(candidate_connect,i);
        for j=1:length(candidate_connect)
            reconstrut_graph(graph_index,1:2)=[i,candidate_connect(j)];
%             reconstrut_graph(graph_index,1:2)=[candidate_connect(j),i];
            graph_index=graph_index+1;
        end
        
%         candidate_connect=[];
%         for j=1:length(temp_r)
%             candidate_connect=[candidate_connect,find(recons_s(j,condidate_index)>0)];
%         end
%         candidate_connect=unique(candidate_connect);
%         for j=1:length(candidate_connect)
%             reconstrut_graph(graph_index,1:2)=[candidate_connect(j),i];
%             graph_index=graph_index+1;
%         end
        
    else
        candidate_connect=[];
        for j=1:length(temp_p)
%             candidate_connect=[candidate_connect,find(recons_s(temp_p(j),:)>0)];
            candidate_connect=[candidate_connect,find(recons_s(temp_p(j),:)<0)];

            
%             temp_revese=find(recons_s(temp_p(j),:)<0);
            temp_revese=find(recons_s(temp_p(j),:)>0);
            temp_revese=setdiff(temp_revese,i);
            temp_revese=temp_revese(find(temp_revese<=length(Ecoli_Model.rxns)));
            temp_revese=temp_revese(find(Ecoli_Model.rev(temp_revese)==1));
            if ~isempty(temp_revese)
%                 candidate_connect=[candidate_connect,temp_revese];
            end
        end
        candidate_connect=unique(candidate_connect);
        candidate_connect=setdiff(candidate_connect,i);
        for j=1:length(candidate_connect)
%             reconstrut_graph(graph_index,1:2)=[candidate_connect(j),i];
            reconstrut_graph(graph_index,1:2)=[i,candidate_connect(j)];
            graph_index=graph_index+1;
        end
        
        candidate_connect=[];
        for j=1:length(temp_r)
%             candidate_connect=[candidate_connect,find(recons_s(temp_r(j),:)>0)];
            candidate_connect=[candidate_connect,find(recons_s(temp_r(j),:)<0)];
            
%              temp_revese=find(recons_s(temp_r(j),:)<0);
            temp_revese=find(recons_s(temp_r(j),:)>0);
            temp_revese=setdiff(temp_revese,i);
            temp_revese=temp_revese(find(temp_revese<=length(Ecoli_Model.rxns)));
            temp_revese=temp_revese(find(Ecoli_Model.rev(temp_revese)==1));
            if ~isempty(temp_revese)
%                 candidate_connect=[candidate_connect,temp_revese];
            end
        
        end
        candidate_connect=unique(candidate_connect);
        candidate_connect=setdiff(candidate_connect,i);
        for j=1:length(candidate_connect)
%              reconstrut_graph(graph_index,1:2)=[candidate_connect(j),i];
           reconstrut_graph(graph_index,1:2)=[i,candidate_connect(j)];
            graph_index=graph_index+1;
        end
    end
end
file_reacgraph=fopen('reaction_group_graph_iMM904.txt','w');
for i=1:size(reconstrut_graph,1)
    fprintf(file_reacgraph,'%d\t%d\n',reconstrut_graph(i,1),reconstrut_graph(i,2));
end
fclose(file_reacgraph);
%% short path length
reconstrut_graph_adj=zeros(size(recons_s,2));
for i=1:size(reconstrut_graph,1)
    reconstrut_graph_adj(reconstrut_graph(i,1),reconstrut_graph(i,2))=1;
end
for i=1:size(reconstrut_graph_adj,1)
    temp=find(reconstrut_graph_adj(i,:)==0);
    if ~isempty(temp)
        reconstrut_graph_adj(i,temp)=inf;
    end
end
% Dist_shortestPath=[];
% for i=1:size(reconstrut_graph_adj,1)
%     temp=simple_dijkstra(reconstrut_graph_adj,i);
%     Dist_shortestPath(i,1:length(temp))=temp;
% end
load -mat Dist_shortestPath
for i=1:size(Dist_shortestPath,1)
    temp=find(abs(Dist_shortestPath(i,:))==inf);
    if ~isempty(temp)
        Dist_shortestPath(i,temp)=50;
    end
end
for i=1:size(Dist_shortestPath,1)
    for j=i:size(Dist_shortestPath,2)
        temp=min(Dist_shortestPath(i,j),Dist_shortestPath(j,i));
        Dist_shortestPath(i,j)=temp;
        Dist_shortestPath(j,i)=temp;
    end
end

spatralCls_Dist_shortestPath=Dist_shortestPath;
spatralCls_index=[];
for i=1:size(spatralCls_Dist_shortestPath,1)
    temp=find(spatralCls_Dist_shortestPath(i,:)~=0);
    spatralCls_Dist_shortestPath(i,temp)=50./spatralCls_Dist_shortestPath(i,temp);
    temp=find(spatralCls_Dist_shortestPath(i,:)==1);
    if length(temp)+1<size(spatralCls_Dist_shortestPath,1)
        spatralCls_index=[spatralCls_index,i];
    end
    temp=find(spatralCls_Dist_shortestPath(i,:)<7);
    spatralCls_Dist_shortestPath(i,temp)=0;
end
spatucl_cluster_candidate=GCSpectralClust2(spatralCls_Dist_shortestPath(spatralCls_index,spatralCls_index),20,20);
spatucl_cluster_candidate_temp=spatucl_cluster_candidate;
Mbst=CNModul(spatucl_cluster_candidate,spatralCls_Dist_shortestPath(spatralCls_index,spatralCls_index));
spatucl_cluster_candidate=spatucl_cluster_candidate(:,Mbst);

cluster_candidate_temp=Dist_shortestPath(spatralCls_index,spatralCls_index);
cluster_candidate_dis=squareform(cluster_candidate_temp);
cluster_candidate_dis(find(cluster_candidate_dis==1))=50;
cluster_candidate=linkage(cluster_candidate_dis,'average');
cluster_candidate_evaluate=cophenet(cluster_candidate,cluster_candidate_dis);
Incon = inconsistent(cluster_candidate);
T=cluster(cluster_candidate,'cutoff',6.6,'criterion','distance');
set(0,'RecursionLimit',2500);
[H,T1,opuperm]=dendrogram(cluster_candidate,0);
 %% some statistic for max cluster
% anabolism num in max 
subCluster_Set={};
 temp=find(T==1);
 index=1;
 maxcluster_index=1;
 while(~isempty(temp))
     subCluster_Set{index}=temp;
     if length(temp)>length(subCluster_Set{maxcluster_index})
         maxcluster_index=index;
     end
     index=index+1;
     temp=find(T==index);
 end
 cluster_subsystem_map={};
 cluster_index=1;
 for i=1:length(subCluster_Set)
     if length(subCluster_Set{i})>10
         cluster_subsystem_map{cluster_index}=spatralCls_index(subCluster_Set{i});
         cluster_index=cluster_index+1;
     end
 end


 group_graph_file=fopen('group_graph_reactionName.txt','w');
 group_reaction_file=fopen('group_reactionName.txt','w');
 for i=1:length(spatralCls_index)
     if spatralCls_index(i)<=length(Ecoli_Model.rxns)
         fprintf(group_graph_file,'%s\t%s\n',Ecoli_Model.rxns{spatralCls_index(i)},Ecoli_Model.grRules{spatralCls_index(i)});
     else
         temp_index=spatralCls_index(i);
         flag=0;
         for j=1:length(group_raction_set)
             temp=group_raction_set{j,8};
             for k=1:length(temp)
                 if temp(k)==temp_index
                     temp_index=j;
                     flag=1;
                     break;
                 end
             end
             if flag==1;
                 break;
             end
         end
         temp=group_raction_set{temp_index,3};
         for j=1:length(temp)-1
            fprintf(group_reaction_file,'%s, ',Ecoli_Model.rxns{temp(j)});
         end
         fprintf(group_reaction_file,'%s\t',Ecoli_Model.rxns{temp(length(temp))});
         fprintf(group_reaction_file,'%s\t',group_raction_set{temp_index,7});
         temp=group_raction_set{temp_index,8};
         for j=1:length(temp)-1
            fprintf(group_reaction_file,'%d, ',temp(j));
         end
         fprintf(group_reaction_file,'%d\n',temp(length(temp)));
     end
 end
 fclose(group_graph_file);
 fclose(group_reaction_file);
 
% enzymeStatistic in max cluster
maxcluster_temp_index=spatralCls_index(find(T==maxcluster_index));
% maxcluster_temp_index=maxcluster_temp_index(find(maxcluster_temp_index<=length(enzyme_reaction_match)));
% enzyme_index_Mcluster=unique(enzyme_reaction_match(maxcluster_temp_index));
% subsystem_map
 cofactor_biosys={'Transport Outer Membrane','Transport Inner Membrane','Transport Outer Membrane Porin'};
 cofactor_biosys_set=[];
 for i=1:length(Ecoli_Model.subSystems)
     for j=1:length(cofactor_biosys)
        if(strcmp(Ecoli_Model.subSystems{i},cofactor_biosys{j}))
            cofactor_biosys_set=[cofactor_biosys_set,i];
            break;
        end
     end
 end
 cofactor_subsystem_map=intersect(maxcluster_temp_index,cofactor_biosys_set);
 maxcluster_temp_index=spatralCls_index(find(T==maxcluster_index));
 cluster_max_index=setdiff(maxcluster_temp_index,cofactor_subsystem_map);
 subnetwork_properity(cluster_max_index)=1;
 enzyme_index_Mcluster=unique(enzyme_reaction_match(cluster_max_index(find(cluster_max_index<=length(enzyme_reaction_match)))));
 
 cluster_max_raction=reconstrut_graph_adj(cluster_max_index,cluster_max_index);
 cluster_max_raction(find(cluster_max_raction==inf))=0;

   
% the rest reaction graph without input 
rest_reaction_index=setdiff(spatralCls_index,maxcluster_temp_index);
rest_input_reaction=intersect(rest_reaction_index,cofactor_biosys_set);
rest_reaction_index=setdiff(rest_reaction_index,rest_input_reaction);
enzyme_index_rest=unique(enzyme_reaction_match(rest_reaction_index(find(rest_reaction_index<=length(enzyme_reaction_match)))));
rest_reaction_set=reconstrut_graph_adj(rest_reaction_index,rest_reaction_index);
rest_reaction_set(find(rest_reaction_set==inf))=0;
input_prephery=[];
 for i=1:size(rest_reaction_set,1)
     if sum(rest_reaction_set(i,:))+sum(rest_reaction_set(:,i))==0
         input_prephery=[input_prephery,rest_reaction_index(i)];
     end
 end
 input_prephery=input_prephery(find(Ecoli_Model.rev(input_prephery)==1));
 input_prephery_reactionName={};
 for i=1:length(input_prephery)
     input_prephery_reactionName=Ecoli_Model.rxns{input_prephery(i)};
 end
 rest_reaction_index=setdiff(rest_reaction_index,input_prephery);
 rest_reaction_set=reconstrut_graph_adj(rest_reaction_index,rest_reaction_index);
 rest_reaction_set(find(rest_reaction_set==inf))=0;
%% the three sub-network model

%% strong compenent decomposition
adjacentGraph=zeros(size(recons_s,2));
for i=1:size(reconstrut_graph,1)
    adjacentGraph(reconstrut_graph(i,1),reconstrut_graph(i,2))=1;
end

 [strong_connect_component_num,strong_connect_component] = graphconncomp(sparse(adjacentGraph));
 topology_matrix=zeros(strong_connect_component_num,strong_connect_component_num);
 strong_connect_component_sort=[];
 topology_rest_attribute=[];
%  for i=1:strong_connect_component_num
%     strong_connect_component_sort=[strong_connect_component_sort,length(temp)];
%  end
%   strong_connect_component_sort=sort(strong_connect_component_sort);
  
 for i=1:strong_connect_component_num
     temp=find(strong_connect_component==i);
     
     if length(temp)>1
%          temp_row=sum(rest_reaction_set(:,temp),2);
         temp_colmn=sum(adjacentGraph(temp,:),1);
%          temp_row_index=find(temp_row~=0);
         temp_colmn_index=find(temp_colmn~=0);
%          topology_matrix(strong_connect_component(temp_row_index),i)=1;
         topology_matrix(i,strong_connect_component(temp_colmn_index))=1;
         topology_rest_attribute=[topology_rest_attribute,4];
     else
        link_edge=find(adjacentGraph(temp,:)==1);
        topology_matrix(i,strong_connect_component(link_edge))=1;
%         link_edge=find(rest_reaction_set(:,temp)==1);
%         topology_matrix(strong_connect_component(link_edge),i)=1;
%         topology_rest_attribute=[topology_rest_attribute,rest_reaction_set_anabolism(temp)];
     end
 end
 for i=1:size(topology_matrix,1)
     topology_matrix(i,i)=0;
 end
 filename=fopen('topologygraph_iMM904.txt','w');
 for i=1:size(topology_matrix,1)
     temp=find(topology_matrix(i,:)~=0);
     if ~isempty(temp)
         for j=1:length(temp)
            fprintf(filename,'%d\t%d\n',i,temp(j));
         end
     end
 end
 fclose(filename);
 order_topology = graphtopoorder(sparse(topology_matrix'));
 
 %% remove reversible edge
 for i=1:size(adjacentGraph,1)
     temp = find(adjacentGraph(i,:)==1);
     if isempty(temp)
         continue;
     end
     for j=1:length(temp)
         adjacentGraph(temp(j),i)=0;
     end
 end
  filename=fopen('norepeatedge_iAF987.txt','w');
 for i=1:size(adjacentGraph,1)
     temp=find(adjacentGraph(i,:)~=0);
     if ~isempty(temp)
         for j=1:length(temp)
            fprintf(filename,'%d\t%d\n',i,temp(j));
         end
     end
 end
 fclose(filename);