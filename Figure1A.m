
clear; clc;
%% Figure 1A; epicenters mapping
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/freesurfer_matlab'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/ENIGMA_matlab'))
addpath(genpath('/media/shuang/data/repository/matlab_packages/plotSurfaceROIBoundary-master'))
addpath(genpath('/media/shuang/data/repository/matlab_packages/brainstat_matlab'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/repository/matlab_packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load cortical thickness matrix and covariates
load('/media/shuang/data/repository/data/CT400.mat'); % cortical thickness
eos = [1:95]; hc = [96:194]; % EOS: early-onset schizophrenia, 95 subjects; HC: healthy controls, 99 subjects
load('/media/shuang/data/repository/data/regressors.mat'); 
%1:95 EOS; 96:end HC --age sex

%% diagnose effect
% ------------ variates -----------------
age=regressors(:,1);
age2=age.^2;
gender = regressors(:,2);
for i =1:size(regressors,1)
    if gender(i)==1
        gender1{i,1}= 'M';
    elseif gender(i)==2
        gender1{i,1}='F';
    end
end
for i=1:95
    diag{i,1}='S';
end
for i=96:length(age)
    diag{i,1}='H';
end
demographics = table(regressors(:,1),age2,diag,gender1,...
    'VariableNames',{'age','age2', 'diag','sex'});
term_sex = FixedEffect(demographics.sex);
term_age = FixedEffect(demographics.age, 'Age');
term_age2 = FixedEffect(demographics.age2, 'Age2');
term_diag = FixedEffect(demographics.diag);
% ------------------ end --------------------------------

model_main = term_sex + term_age + term_diag;

% %--------------validation analysis-------------
% model_main = term_sex  + term_age + term_age2 + term_diag; 
% %-----------------end---------------------------
 
contrast_diag = (demographics.diag == "S")-(demographics.diag == "H");
slm_diag = SLM( ...
    model_main, ...
    contrast_diag,...
    'cluster_threshold', 0.05, ...
    'correction',{'fdr'});
slm_diag.fit(CT400');
p=tcdf(-abs(slm_diag.t),slm_diag.df);
FDRvalue = mafdr(p,'BHFDR',true); % 
%-----------------show diagnose effect-----------------
obj=plot_hemispheres([(slm_diag.t)', (slm_diag.t.*(FDRvalue<0.05))'], ...
             {surf_lh,surf_rh}, 'parcellation', labeling, ...
             'labeltext',{'diag t', 'FDRDt'});
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-4, 4]) 
%------------------end----------------------------------

tmap=slm_diag.t';

%%---- save CT tmap for diagnose effect -----
mask=ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
fdr05=(FDRvalue<0.05)';
for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
    mask.cdata(find(mask.cdata==i))=fdr05(i);
end
ciftisavereset(mask,'/media/shuang/data/repository/result/CTdiagFDR.dlabel.nii',path_wb_command);

mask=ciftiopen('/media/shuang/data/repository/templates/SCON003_400.pscalar.nii',path_wb_command);
mask.cdata=tmap;
ciftisavereset(mask,'/media/shuang/data/repository/result/CTdiag.pscalar.nii',path_wb_command);
%%------end-------------

%% Functional epicenters mapping
load('/media/shuang/data/repository/data/FC_connectome.mat'); % FC matrix
% negative to 0;z-transform 
for i = 1:177
    zfc(:,:,i)=reshape(normalize(reshape(connectome(:,:,i),400*400,1),'range',[-1,1]),400,400); % normalize z to 0 1 range
end
zfc(find(zfc<=0))=0;
for i=1:177
    zfc(:,:,i)=atanh(zfc(:,:,i)); % Fisher-z transformed
end
FC=mean(zfc(:,:,87:end),3); % HC---normative functional connectome
FC(isinf(FC))=1;

%--------------plot heatmap of FC matrix---------------
figure; heatmap(FC,'ColorLimits',[0 0.6],'GridVisible','off'); 
colormap(brewermap([],"Greys"))

%--------------plot FC degree centrality--------------
obj=plot_hemispheres(sum(FC,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','FC degree');
obj.colormaps(brewermap([],"Greys"))
obj.colorlimits([10, 80]) 
%----------------- end -----------------------------

%----------------calculating functional epicenters---------
for i=1:400
    t1=tmap;t1(i)=[];
    t2=FC(:,i);t2(i)=[];
    epifc(i)=corr(t1,t2);
end

obj=plot_hemispheres(epifc',{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','Functional epicenters');
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.2, 0.2])

% ------------------ spin test ---------------------
load('/media/shuang/data/repository/data/perm_id.mat');
spins = perm_id;          % spatial autocorrelation-preserving permutation assignments
nspins = 10000;                 % number of permutations ("spins")
epifc_perms = zeros(400,nspins);
for k=1:nspins
    for i=1:400
        epifc_perms(i,k)=corr(tmap,FC(spins(:,k),i));
    end
end

for i = 1:400 % one-tailed
    if epifc(i)>0
        spinp(i) = (1+(nnz(find(epifc_perms(i,:)>=epifc(i)))))/(1+nspins);
    else
        spinp(i) = (1+(nnz(find(epifc_perms(i,:)<=epifc(i)))))/(1+nspins);
    end
end

spinp05=spinp<0.05;
obj=plot_hemispheres([epifc',epifc'.*spinp05],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'epicenterFC','spin0.05'});
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.2, 0.2])

mask=ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
FCspin05=spinp05;
for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
    mask.cdata(find(mask.cdata==i))=FCspin05(i);
end
ciftisavereset(mask,'/media/shuang/data/repository/result/EpiFCspin05.dlabel.nii',path_wb_command);

mask=ciftiopen('/media/shuang/data/repository/templates/SCON003_400.pscalar.nii',path_wb_command);
mask.cdata=epifc';
ciftisavereset(mask,'/media/shuang/data/repository/result/epifc.pscalar.nii',path_wb_command);

%% structural connectome
load('/media/shuang/data/repository/data/SC_connectome.mat');
SC=mean(Sconnectome,3); %HC
SC=SC+eye(400);

%------ plot normative structural connectome
figure; heatmap(SC,'ColorLimits',[0 0.6],'GridVisible','off'); 
colormap(brewermap([],"Greys"))

obj=plot_hemispheres(sum(SC,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','SC degree');
obj.colormaps(brewermap([],"Greys"))
obj.colorlimits([1, 40]) 
%------ end -------


for i=1:400
    t1=tmap;t1(i)=[];
    t2=SC(:,i);t2(i)=[];
    episc(i)=corr(t1,t2);
end

obj=plot_hemispheres(episc',{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','SCepic');
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.2, 0.2])

%-------------spin test-------------
episc_perms = zeros(400,nspins,3);
for k=1:nspins
    for i=1:400
        episc_perms(i,k)=corr(tmap,SC(spins(:,k),i));
    end
end

for i = 1:400 % one-tailed
    if episc(i)>0
        spinp(i) = (1+(nnz(find(episc_perms(i,:)>=episc(i)))))/(1+nspins);
    else
        spinp(i) = (1+(nnz(find(episc_perms(i,:)<=episc(i)))))/(1+nspins);
    end
end
spinp05=spinp<0.05;

obj=plot_hemispheres([episc',episc'.*spinp05],{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'DiagSC','spin0.05 Inter'});
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.2, 0.2])

mask=ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
SCspin05=spinp05;
for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
    mask.cdata(find(mask.cdata==i))=SCspin05(i);
end
ciftisavereset(mask,'/media/shuang/data/repository/result/EpiSCspin05.dlabel.nii',path_wb_command);

mask=ciftiopen('/media/shuang/data/repository/templates/SCON003_400.pscalar.nii',path_wb_command);
mask.cdata=episc';
ciftisavereset(mask,'/media/shuang/data/repository/result/episc.pscalar.nii',path_wb_command);

