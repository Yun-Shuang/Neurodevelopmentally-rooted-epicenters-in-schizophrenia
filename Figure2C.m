
clear;clc;
%% Figure 2C; epicenters cognitive relevance---Neurosynth decoding

addpath(genpath('/media/shuang/data/repository/matlab_packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/freesurfer_matlab'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/repository/matlab_packages/ENIGMA_matlab'))
addpath(genpath('/media/shuang/data/repository/matlab_packages/meta'))

%% load data
load('/media/shuang/data/repository/public_data/CT_psych_gradients.mat');
load('/media/shuang/data/repository/result/epifc.mat');
load('/media/shuang/data/repository/result/episc.mat');
epifc_surf = parcel_to_surface(epifc,'schaefer_400_conte69');
fcepi = surface_to_parcel(epifc_surf, 'schaefer_400_conte69');
episc_surf = parcel_to_surface(episc,'schaefer_400_conte69');
scepi = surface_to_parcel(episc_surf, 'schaefer_400_conte69');

cd /media/shuang/data/repository/public_data/neurosynth_z_values;
BH =     {'action_association-test_z_lh.shape.gii'
    'action_association-test_z_rh.shape.gii'
    'affective_association-test_z_lh.shape.gii'
    'affective_association-test_z_rh.shape.gii'
    'attention_association-test_z_lh.shape.gii'
    'attention_association-test_z_rh.shape.gii'
    'auditory_association-test_z_lh.shape.gii'
    'auditory_association-test_z_rh.shape.gii'
    'autobiographical_memory_association-test_z_lh.shape.gii'
    'autobiographical_memory_association-test_z_rh.shape.gii'
    'cognitive_control_association-test_z_lh.shape.gii'
    'cognitive_control_association-test_z_rh.shape.gii'
    'emotion_association-test_z_lh.shape.gii'
    'emotion_association-test_z_rh.shape.gii'
    'episodic_memory_association-test_z_lh.shape.gii'
    'episodic_memory_association-test_z_rh.shape.gii'
    'eye_movement_association-test_z_lh.shape.gii'
    'eye_movement_association-test_z_rh.shape.gii'
    'face_association-test_z_lh.shape.gii'
    'face_association-test_z_rh.shape.gii'
    'inhibition_association-test_z_lh.shape.gii'
    'inhibition_association-test_z_rh.shape.gii'
    'language_association-test_z_lh.shape.gii'
    'language_association-test_z_rh.shape.gii'
    'motor_association-test_z_lh.shape.gii'
    'motor_association-test_z_rh.shape.gii'
    'multisensory_association-test_z_lh.shape.gii'
    'multisensory_association-test_z_rh.shape.gii'
    'pain_association-test_z_lh.shape.gii'
    'pain_association-test_z_rh.shape.gii'
    'reading_association-test_z_lh.shape.gii'
    'reading_association-test_z_rh.shape.gii'
    'reward_association-test_z_lh.shape.gii'
    'reward_association-test_z_rh.shape.gii'
    'semantics_association-test_z_lh.shape.gii'
    'semantics_association-test_z_rh.shape.gii'
    'social_cognition_association-test_z_lh.shape.gii'
    'social_cognition_association-test_z_rh.shape.gii'
    'verbal_association-test_z_lh.shape.gii'
    'verbal_association-test_z_rh.shape.gii'
    'visual_association-test_z_lh.shape.gii'
    'visual_association-test_z_rh.shape.gii'
    'visual_perception_association-test_z_lh.shape.gii'
    'visual_perception_association-test_z_rh.shape.gii'
    'visuospatial_association-test_z_lh.shape.gii'
    'visuospatial_association-test_z_rh.shape.gii'
    'working_memory_association-test_z_lh.shape.gii'
    'working_memory_association-test_z_rh.shape.gii'};
LH_ind = find(endsWith(BH, 'lh.shape.gii'));
LH = BH(LH_ind);
RH_ind = find(endsWith(BH, 'rh.shape.gii'));
RH = BH(RH_ind);

%% bin neurosynth data
load('/media/shuang/data/repository/templates/ENIGMA_parcellations/schaefer_400_conte69.csv');
nb=40;
tmp = nan(24,nb);
for j = 1:24
    L = gifti(strcat('/media/shuang/data/repository/public_data/neurosynth_z_values/', LH{j}));
    R = gifti(strcat('/media/shuang/data/repository/public_data/neurosynth_z_values/', RH{j}));
    funcmap= [L.cdata;R.cdata];
    bins1   = quantileranks(fcepi,nb); % binned into 1:40 bins and allocated to 400 parcels
    epi_fc = zeros(1,64984);
    for i = 1:length(bins1)
        a = i-1; 
        epi_fc(:,schaefer_400_conte69==a) = bins1(i);
    end
    for i = 1:nb
        tmp(j,i) = mean(funcmap(epi_fc==i)); %size: 24 x 40
    end
end

%zscore binned epicenter map
epifc_binned = nan(1,nb);
for i = 1:nb
    epifc_binned(i) = mean(fcepi(bins1==i));
end
tmpz2= zscore(tmp,[],2);
tmpz22 = tmpz2;
tmpz2(tmpz2<0.5) = nan;
weighted_mean = nanmean(tmpz2.*(epifc_binned),2);
[~,b,c] = unique(weighted_mean);

tmp2 = nan(24,nb);
epi_binned = nan(1,nb);
for j = 1:24
    L = gifti(strcat('/media/shuang/data/repository/public_data/neurosynth_z_values/', LH{j}));
    R = gifti(strcat('/media/shuang/data/repository/public_data/neurosynth_z_values/', RH{j}));
    funcmap= [L.cdata;R.cdata];
    bins1   = quantileranks(scepi,nb);
    epi_fc = zeros(1,64984);
    for i = 1:length(bins1)
        a = i-1;
        epi_fc(:,find(schaefer_400_conte69==a)) = bins1(i);
    end
    for i = 1:nb
        tmp2(j,i) = mean(funcmap(find(epi_fc==i)));
    end
end
for i = 1:nb
    epi_binned(i) = mean(scepi(bins1==i));
end
tmpz= zscore(tmp2,[],2);
tmpzz = tmpz;
tmpz(tmpz<0.5) = nan;
weighted_mean = nanmean(tmpz.*(epi_binned),2);
[TT,d,e] = unique(weighted_mean);

j=0;
for i = 1:24
    j=j+1;
    words{j} = LH{j};
    worrs3{j} = words{j}(1:end-32);
    words_final{j} = strrep(worrs3{j},'_',' ');
    wrd{j} = num2str(j);
end
eucd= sqrt(c.^2+d.^2);

%%  color palette
epifc_plot = epifc;
episc_plot = episc;
n=400;
g = load('/media/shuang/data/repository/matlab_packages/meta/cb/cbrewer/colorbrewer.mat');
colorbrewer = g.colorbrewer;

epifc_peak = max(epifc_plot);
episc_min = min(episc_plot);
episc_max = max(episc_plot);

colourness = nan(n,3);
for ii = 1:length(epifc_plot)
    
    if epifc_plot(ii) > 0
        colourness(ii,2) = 1 - (epifc_peak - epifc_plot(ii));
    else
        colourness(ii,2) = 1 - (epifc_peak + abs(epifc_plot(ii)));
    end
    if episc_plot(ii) > 0
        colourness(ii,1) = 1 - (episc_max- episc_plot(ii));
        colourness(ii,3) = 1 - (abs(episc_min) + episc_plot(ii));
    else
        colourness(ii,1) = 1 - (episc_max + abs(episc_plot(ii)));
        colourness(ii,3) = 1 - (abs(episc_min) - abs(episc_plot(ii)));
    end
    
end
% rescale colorbar 
colourness_rescale = nan(n,3);
colourness_vertices = nan(64984,3);
for col = 1:3
    colourness_rescale(:,col) = rescale(colourness(:,col), 0, 1);
    colourness_vertices(:,col)=parcel_to_surface(colourness_rescale(:,col),'schaefer_400_conte69');
end

%compute densities for density plot next to scatters
[x, y] = ksdensity(epifc_plot);
y2(1,:) = x; x2(1,:) = y;
[y2(2,:), x2(2,:)] = ksdensity(episc_plot);

%% scatter
f= figure;
a(1) = axes('position', [0 0.1 0.7 0.7]);
scatter(epifc_plot, episc_plot, 20, colourness_rescale,'filled');
axis([min(epifc_plot)-0.01 max(epifc_plot)+0.01 min(episc_plot)-0.01 max(episc_plot)+0.01])
set(gca,'xtick',[],'ytick',[])
axis square

a(2) = axes('position', [0.65 0.1 0.07 0.7]); %add densities
patch(y2(2,:),x2(2,:), ones(1,length(y2(2,:)))); axis off;
colormap(a(2), [190, 190, 190]/250)
ylim([min(episc_plot) max(episc_plot)])

a(3) = axes('position', [0.1 0.85 0.5 0.07]);
patch(x2(1,:), y2(1,:), ones(1,length(y2(1,:)))); axis off;
colormap(a(3), [190, 190, 190]/250)
xlim([min(epifc_plot) max(epifc_plot)])

hold on
a(4) = axes('position', [0 0.1 0.7 0.7]); col = [0,0,0]; %add cognitive function terms
x = c; y = e; scatter(x,y,40,col,'filled');
txty = words_final;
dx = 0.85;
dy = 0.2; 
text(x+dx, y+dy, txty,'FontSize',11, 'FontName', 'Aril');
axis square
set(gca,'XTick',[],'YTick',[],'color','none');xlabel('fcepi');ylabel('scepi')
hold off
