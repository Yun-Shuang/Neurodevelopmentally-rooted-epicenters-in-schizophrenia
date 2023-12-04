
clear;clc;
%% Figure3: PLS analyses between gene expression maps and disease epicenters
addpath(genpath('/media/shuang/data/repository/matlab_packages/brainstat_matlab'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/hansen_genescognition-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrewerMap-master'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/repository/matlab_packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% Figure 3A: PLS analysis
load('/media/shuang/data/repository/public_data/expression.mat');
X=zscore(expression);
load('/media/shuang/data/repository/result/epifc.mat');
load('/media/shuang/data/repository/result/episc.mat');
Y=zscore([epifc, episc]);
figure; heatmap(X,'ColorLimits',[-4 4],'GridVisible','off'); colormap(flipud(brewermap([],"PiYG")))
figure; heatmap(Y,'ColorLimits',[-3 3],'GridVisible','off'); colormap(flipud(brewermap([],"PiYG")))

nnodes = size(X,1);
nterms = size(Y,2);
ngenes = size(X,2);
% behav pls
option.method = 3;
option.num_boot = 10000;
option.num_perm = 0;              
option.stacked_behavdata = Y;
exp{1} = X;
result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
load('/media/shuang/data/repository/result/PLSresult.mat');

%% spin test: find significant latent variable
load('/media/shuang/data/repository/data/perm_id.mat');
spins = perm_id;          % spatial autocorrelation-preserving permutation assignments
nspins = 10000;                 % number of permutations ("spins")
s_spins = zeros(nterms,nspins); % singular values
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
for k = 1:nspins    
    option.stacked_behavdata = Y(spins(:,k),:);  % permute neurosynth matrix
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0); 
    [r,c] = size(datamatsvd);
    if r <= c
        [pu, sperm, pv] = svd(datamatsvd',0);
    else
        [pv, sperm, pu] = svd(datamatsvd,0);
    end
    %  rotate pv to align with the original v
    rotatemat = rri_bootprocrust(result.v,pv);
    %  rescale the vectors
    pv = pv * sperm * rotatemat;
    sperm = sqrt(sum(pv.^2));
    s_spins(:,k) = sperm;
end
sprob = zeros(nterms,1); % p-value for each latent variable
for k = 1:nterms % get permuted (via spin test) p-values
    sprob(k) = (1+(nnz(find(s_spins(k,:)>=result.s(k)))))/(1+nspins); % terms
end  

%% plot effect size for PLS
ES=(result.s).^2; ES=ES./sum(ES);
ES_null=(s_spins).^2;
for i=1:size(ES_null,2)
    ES_null(:,i)=ES_null(:,i)./sum(ES_null(:,i));
end
ESsprob = zeros(nterms,1); % p-value for each latent variable
for k = 1:nterms % get permuted (via spin test) p-values
    ESsprob(k) = (1+(nnz(find(ES_null(k,:)>=ES(k)))))/(1+nspins); % terms
end

figure;
scatter(1:2,ES,30,'filled')
hold on
boxplot(squeeze(ES_null)')
set(gca,'xticklabel')
xtickangle(90)
ylim([0 1])
title('pls perm')

%% cross validation the covarying patterns between genes and epicenters
lv=1; % use the first variable for only it is significant
figure; plot(result.usc(:,lv),result.vsc(:,lv),'k.','markersize',20);
ylim([-3 3]);xlim([-125 125]);
hold on; plot(result.usc(:,lv),polyval(polyfit(result.usc(:,lv),result.vsc(:,lv),1),result.usc(:,lv)),'r')

exp{1} = X;
[rtrain,rtest] = fcn_crossval_pls_brain_obvs(exp,Y,100,0.75,lv,coords); 
nperm = 10000;
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
null_corr=zeros(nperm,1);
for m = 1:nperm % for each permutation
    option.stacked_behavdata = Y(perm_id(:,m),:); % rows of neurosynth matrix has been permuted while preserving spatial autocorrelation
    result = pls_analysis(exp, nnodes, 1, option);              % this is the null PLS result
    null_corr(m)=corr(result.usc(:,1),result.vsc(:,1));
end
COsprob = zeros(nterms,1); % p-value for each latent variable
for k = 1:nterms % get permuted (via spin test) p-values
    COsprob(k) = (1+(nnz(find(null_corr(k,:)>=rtest(k)))))/(1+nperm); % terms
end

figure;
boxplot(null_corr)
hold on; 
boxplot([rtest rtrain])
ylim([0 0.8])
set(gca,'xticklabel',{'test','train'})
ylabel('correlation')
title('cross validation')

%% Figure 3B: plot gene and epicenter scores
obj=plot_hemispheres(result.usc(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','Gene scores');
obj.colormaps([0.7 0.7 0.7;flipud(brewermap([],'PiYG'))])
obj.colorlimits([-125, 125]) % Gene scores for PLS1

obj=plot_hemispheres(result.vsc(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','Epicenter scores');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],'PiYG'))])
obj.colorlimits([-3, 3]) % epicenter scores

%%  Figure 3C: get significant genes and their loadings
load('/media/shuang/data/repository/public_data/gene_names.mat'); % all gene names
% spin test to generate null distributions of gene loadings
option.method = 3;
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
nperm = 10000;
gload_null = zeros(ngenes,nperm);
for m = 1:nperm % for each permutation
    clear option.stacked_behavdata
    option.stacked_behavdata = Y(perm_id(:,m),:); 
    result = pls_analysis(exp, nnodes, 1, option);
    for k = 1:ngenes
        gload_null(k,m) = corr(expression(:,k),result.usc(:,1)); % null loadings
    end
end
sprob = zeros(size(gload_null,1),1); 
for k = 1:size(gload_null,1) % get permuted (via spin test) p-values
    if gload(k)>0
        sprob(k) = (1+(nnz(find(gload_null(k,:)>=gload(k)))))/(1+nperm);
    else
        sprob(k) = (1+(nnz(find(gload_null(k,:)<=gload(k)))))/(1+nperm);
    end
end  
genesid=sprob<0.05;

clear result
load('/media/shuang/data/repository/result/PLSresult.mat');
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expression(:,k),result.usc(:,1));
end
gload=gload.*genesid;
ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos,'descend'); % sorted
[~,Ineg] = sort(gload_neg); % sorted 
threshold = 1; % all genes: could set it lower to show only top genes
gpos_idx = Ipos(1:floor(threshold*length(gload_pos)));
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));

pname=gene_name(ipos(gpos_idx))';pp = cell2table(cellstr(pname));
pp.('weights')=abs(gload(ipos(gpos_idx)));
writetable(pp,'/media/shuang/data/repository/result/PgenesymbSIG.csv'); % used for metascape input
clear pname pp
pname=gene_name(ineg(gneg_idx))';
pp = cell2table(cellstr(pname));
pp.('weights')=abs(gload(ineg(gneg_idx))); %weights is gene loadings
writetable(pp,'/media/shuang/data/repository/result/NgenesymbSIG.csv');

%% cell types
load('/media/shuang/data/repository/public_data/cell_types.mat');
g=[1:15631]';
% find empirical cell type ratio
ctd_ratios = zeros(ntypes,2);
for k = 1:length(ctd_ratios)
    ctd_ratios(k,1) = length(intersect(g(ipos(gpos_idx)),cell2mat(ctd_idx{k})))/length(gpos_idx);
    ctd_ratios(k,2) = length(intersect(g(ineg(gneg_idx)),cell2mat(ctd_idx{k})))/length(gneg_idx);
end

% null model from random gene set
n = 10000;
ctd_null = zeros(ntypes,2,n);
for k = 1:n 
    y = datasample([1:length(g)],length(gpos_idx),'Replace',false);                      
    for j = 1:ntypes                                                                     
        ctd_null(j,1,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gpos_idx); 
    end
    % negative nulls
    y = datasample([1:length(g)],length(gneg_idx),'Replace',false);                      
    for j = 1:ntypes                                                                     
        ctd_null(j,2,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gneg_idx); 
    end
end
% get p-values with two-tailed significance test
p_ctd = ctd_ratios - mean(ctd_null,3); % mean centre
p_null = ctd_null - mean(ctd_null,3);
pval = zeros(ntypes,2);
for k = 1:ntypes
    pval(k,1) = (1+(nnz(find(abs(p_null(k,1,:)) >= abs(p_ctd(k,1))))))/(n+1); % pval for positive gene set
    pval(k,2) = (1+(nnz(find(abs(p_null(k,2,:)) >= abs(p_ctd(k,2))))))/(n+1); % pval for negative gene set
end
for k = 1:2
    pval(:,k) = mafdr(pval(:,k),'BHFDR',true); % FDR correction
end

% plot empirical ratios as points and null models as boxplots
o = [1,3,6,7,2,4,5]; % order in which cell types appear
op=pval(o,:);
figure;
subplot(1,2,1) % specific cell type expression for positive gene set
scatter(1:ntypes,ctd_ratios(o,1),30,'filled')
hold on
boxplot(squeeze(ctd_null(o,1,:))')
set(gca,'xticklabel',C(o))
xtickangle(90)
ylim([0 0.17])
title('pos')

subplot(1,2,2) % specific cell type expression for negative gene set
scatter(1:ntypes,ctd_ratios(o,2),30,'filled')
hold on
boxplot(squeeze(ctd_null(o,2,:))')
set(gca,'xticklabel',C(o))
xtickangle(90)
ylim([0 0.17])
title('neg')    
