
clear; clc;
%% Figure 4: Disorders and evolution-associated analyses
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/brainstat_matlab'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/repository/matlab_packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load significant genes
load('/media/shuang/data/repository/public_data/gene_names.mat'); % all gene names
celltypes = table2cell(readtable('/media/shuang/data/repository/result/PgenesymbSIG.csv')); % load specific cell type expression
Pgenes = cellstr(celltypes(:,1)); Pw=cell2mat(celltypes(:,2));
celltypes = table2cell(readtable('/media/shuang/data/repository/result/NgenesymbSIG.csv')); % load specific cell type expression
Ngenes = cellstr(celltypes(:,1));Nw=cell2mat(celltypes(:,2));
sigg=zeros(15631,2);
for i=1:size(Pgenes,1)
    sigg(find(strcmp(gene_name,Pgenes(i))),1)=Pw(i);
end
for i=1:size(Ngenes,1)
    sigg(find(strcmp(gene_name,Ngenes(i))),2)=Nw(i);
end

%% Figure 4A: find overlap genes with AHBA
% 6 disorders: ASD; SZ; BD; MDD; Alcoholism; IBD
T = table2cell(readtable('/media/shuang/data/repository/public_data/DGE_6disorders/genesDGEfdr.csv')); % load entrezID of genes
overlapgene = cell(size(gene_name));
idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name(k), T(:,1))                                      
        overlapgene{k} = cell2mat(T(find(strcmp(gene_name(k), T(:,1))),2:end));  
        idx = [idx;k];                                                     % also store the index of the gene
    end
end
overlapgene=overlapgene';
dge=zeros(15631,13);
rep=[];
for i=1:size(overlapgene,1)
    temp=cell2mat(overlapgene(i));
    if isempty(temp)
        dge(i,:)=zeros(1,13);
    elseif size(temp,1)~=1
        temp(isnan(temp))=0;
        dge(i,1)=temp(1,1);
        for j=2:size(temp,2) % some genes might have repeated dge values
            if temp(1,j)&&temp(2,j)
                dge(i,j)=mean(temp(:,j));
            else
                dge(i,j)=temp(1,j)+temp(2,j);
            end
        end
        rep=[rep;i];
    else
        dge(i,:)=temp;
    end
end

%% correlated gene loadings with DGE values
pdge=[sigg(logical(sigg(:,1)),1),dge(logical(sigg(:,1)),2:end)]; % 298 gene weights and dge values
ndge=[-sigg(logical(sigg(:,2)),2),dge(logical(sigg(:,2)),2:end)]; % 493
tem=[pdge;ndge];
nam=[gene_name(logical(sigg(:,1)))';gene_name(logical(sigg(:,2)))'];
nam(any(tem==0,2),:)=[];
pdge(any(pdge==0,2),:)=[]; % 272 discard those without dge
ndge(any(ndge==0,2),:)=[]; % 438
gene=[pdge;ndge];
loop=10000; % permutation tests
for i=1:6 % 6 disorders
    temp=find(gene(:,2*i+1)<0.05);
    ind=zeros(size(gene,1),1);ind(temp,1)=1;ind=logical(ind);
    num(i)=sum(ind); % dge p<0.05
    if num(i)<5
        rg(i)=nan;
    else
    [rg(i),p(i)]=corr(gene(ind,1),gene(ind,2*i),'type','Spearman');    
    end
    x=gene(ind,1);y=gene(ind,2*i); % do permutations    
    hue=zeros(length(x),1);
    hue((sign(x)==1)&(sign(y)==1))=1; hue((sign(x)==1)&(sign(y)==-1))=1;
    hue((sign(x)==-1)&(sign(y)==-1))=2; hue((sign(x)==-1)&(sign(y)==1))=2;
    T=table(x,y,hue);
    writetable(T,['/media/shuang/data/repository/result/DGE_6disorders/disorder', num2str(i),'.csv']);
    for j=1:loop
        wei=x(randperm(length(x)));
        perm_r(j,i)=corr(wei,y);
    end
    if rg(i)>0 
        perm_p(i) = (1+(nnz(find(perm_r(:,i)>=rg(i)))))/(1+loop);
    elseif rg(i)<0
        perm_p(i) = (1+(nnz(find(perm_r(:,i)<=rg(i)))))/(1+loop);
    end
end
FDR = mafdr(perm_p,'BHFDR',true); % FDR correction

%% Figure 4B: human accelerated region genes
load('/media/shuang/data/repository/result/geneloadings.mat');
[B,I]=sort(gload(idx),'descend'); % from PLS+ to PLS-
T=table2cell(readtable('/media/shuang/data/repository/public_data/HARs_genes/HARgene.csv'));
overlapgene = cell(size(gene_name));
idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name(k), T(:,1))                                      
        overlapgene{k} = cell2mat(T(find(strcmp(gene_name(k), T(:,1)))));  
        idx = [idx;k];                                                     % also store the index of the gene
    end
end

load('/media/shuang/data/p_02658/Genetic/expression.mat');
load('/media/shuang/data/repository/result/epifc.mat');
load('/media/shuang/data/repository/result/episc.mat');
X=zscore(expression);
overlap=X(:,idx); % 1372 genes
clear r p
for i=1:size(overlap,2)
    [r(i,1),p(i,1)]=corr(epifc,overlap(:,i));
    [r(i,2),p(i,2)]=corr(episc,overlap(:,i));
end

% ------ spin test --------
load('/media/shuang/data/repository/data/perm_id.mat');
spins = perm_id;         
nspins = 10000;
epifc_perms = zeros(size(overlap,2),nspins);
for k=1:nspins
    for i=1:size(overlap,2)
        epifc_perms(i,k)=corr(epifc,overlap(spins(:,k),i));
        episc_perms(i,k)=corr(episc,overlap(spins(:,k),i));
    end
end
for i = 1:size(overlap,2)
    if r(i,1)>0
        spinp(i,1) = (1+(nnz(find(epifc_perms(i,:)>=r(i,1)))))/(1+nspins);
    else
        spinp(i,1) = (1+(nnz(find(epifc_perms(i,:)<=r(i,1)))))/(1+nspins);
    end
end
for i = 1:size(overlap,2) 
    if r(i,2)>0
        spinp(i,2) = (1+(nnz(find(episc_perms(i,:)>=r(i,2)))))/(1+nspins);
    else
        spinp(i,2) = (1+(nnz(find(episc_perms(i,:)<=r(i,2)))))/(1+nspins);
    end
end
spinp05=spinp<0.05;
ttable=gene_name(idx)'; % name fc spinp0.05 sc spinp0.05
ttable(:,2:3)=[r(:,1),spinp05(:,1)];ttable(:,4:5)=[r(:,2),spinp05(:,2)];
ttable=ttable(I,:);
pp = cell2table(cellstr(ttable));
writetable(pp,'/media/shuang/data/repository/result/HARgenesSIG.csv');

obj=plot_hemispheres(X(:,idx(find(r(:,1)==max(r(:,1))))),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'FKBP5 FChrfPLS+'});
obj.colormaps([0.7 0.7 0.7;flipud(brewermap([],'PiYG'))])
obj.colorlimits([-2.4, 2.4])

obj=plot_hemispheres(X(:,idx(find(r(:,1)==min(r(:,1))))),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'CA10 FChrfPLS-'});
obj.colormaps([0.7 0.7 0.7;flipud(brewermap([],'PiYG'))])
obj.colorlimits([-2.4, 2.4])
