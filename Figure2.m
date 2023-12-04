
clear; clc;
%% Figure 2A; Individualized epicenters mapping
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/brainstat_matlab'));
addpath(genpath('/media/shuang/data/software/ENIGMA/matlab/'))
addpath(genpath('/media/shuang/data/software/plotSurfaceROIBoundary-master/'))
addpath(genpath('/usr/local/freesurfer/matlab'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/repository/matlab_packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load cortical thickness matrix and covariates
load('/media/shuang/data/repository/data/CT400.mat'); % cortical thickness
eos = [1:95]; hc = [96:194]; % EOS: early-onset schizophrenia, 95 subjects; HC: healthy controls, 99 subjects
HC = CT400(:,96:end); EOS = CT400(:,1:95);
load('/media/shuang/data/repository/data/regressors.mat'); 
%1:95 EOS; 96:end HC --age sex

%% Deviations of cortical thickness in EOS patients
w=zeros(size(EOS));
for i=1:400
    w(i,:)=(EOS(i,:)-mean(HC(i,:)))./std(HC(i,:));
end

%% Individulized epicenters mapping
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

load('/media/shuang/data/repository/data/SC_connectome.mat');
SC=mean(Sconnectome,3); %HC
SC=SC+eye(400);

for i=1:400
    for j=1:size(w,2) %95 eos
        epifc(i,j)=corr(w(:,j),FC(:,i));
    end
end
for i=1:400
    for j=1:size(w,2)
        episc(i,j)=corr(w(:,j),SC(:,i));
    end
end

%% Figure 2B; Subgroup epicenters mapping
clear epifc episc;

regressors_ori=regressors;
age=regressors_ori(:,1);
group1=age<14; 
group21=(age>=14)&(age<16); group22=age>=16;
child=regressors(group1,:);
early=regressors(group21,:);
late=regressors(group22,:);

% ------------ variates -----------------
group=group1; % childhood
% group=group21; % early adolescence
% group=group22; % late adolescence
clear regressors gender gender1 age diag demographics
HC=CT400(:,95+find(group(hc)));
EOS=CT400(:,find(group(eos)));
regressors=regressors_ori(group,:);
gender = regressors(:,2);
for i =1:size(regressors,1)
    if gender(i)==1
        gender1{i,1}= 'M';
    elseif gender(i)==2
        gender1{i,1}='F';
    end
end
for i=1:size(EOS,2)
    diag{i,1}='S';
end
for i=(size(EOS,2)+1):length(gender)
    diag{i,1}='H';
end
demographics = table(diag,gender1,...
    'VariableNames',{'diag','sex'});
term_sex = FixedEffect(demographics.sex);
term_diag = FixedEffect(demographics.diag);
% ------------------ end --------------------------------

model_main = term_sex + term_diag + term_diag;
%%-------------diag main effect-------------------------
clear slm_diag p mainT
contrast_diag = (demographics.diag == "S")-(demographics.diag == "H");
slm_diag = SLM( ...
    model_main, ...
    contrast_diag,...
    'cluster_threshold', 0.05, ...
    'correction',{'fdr'});
slm_diag.fit(CT400(:,group)');
p=tcdf(-abs(slm_diag.t),slm_diag.df); %
FDRvalue = mafdr(p,'BHFDR',true); 
tmap(:,1)=slm_diag.t'; % childhood
% tmap(:,2)=slm_diag.t'; % early adolescence
% tmap(:,3)=slm_diag.t'; % late adolescence

%% epicenters mapping --functional epicenters
load('/media/shuang/data/repository/result/subageTMAP.mat');
load('/media/shuang/data/repository/result/fc_subage.mat');
load('/media/shuang/data/repository/result/sc_subage.mat');
for i=1:400
    for j=1:size(tmap,2)
        epifc(i,j)=corr(tmap(:,j),fc(:,i,j));
    end
end

obj=plot_hemispheres(epifc(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','Childepicenters FC');
obj.colorlimits([-0.3, 0.3])
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])

obj=plot_hemispheres(epifc(:,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EarlyAdolepicenters FC');
obj.colorlimits([-0.3, 0.3])
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])

obj=plot_hemispheres(epifc(:,3),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','LateAdolepicenters FC');
obj.colorlimits([-0.3, 0.3])
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])

%% structural epicenters
for i=1:400
    for j=1:size(tmap,2)
        episc(i,j)=corr(tmap(:,j),sc(:,i,j));
    end
end

obj=plot_hemispheres(episc(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','SCepic child');
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))]) % for making no hole
obj.colorlimits([-0.3, 0.3])

obj=plot_hemispheres(episc(:,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','SCepic Earlyadol');
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.3, 0.3])

obj=plot_hemispheres(episc(:,3),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','SCepic Lateyadol');
obj.colormaps([0.7,0.7,0.7;flipud(brewermap([],"RdBu"))])
obj.colorlimits([-0.3, 0.3])

