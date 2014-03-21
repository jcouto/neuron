clear all
close all
clc
% if REMOVE_OUTLIER_TRIALS is selected, the trials that are +- one standard
% deviation of the statistics of the overall firing rates are removed.
REMOVE_OUTLIER_TRIALS = 1;
% Extract all spiketimes
folders = dir();
folders = {folders([folders.isdir]).name};
folders([1,2]) = [];
trials = cell(1,length(folders));
t0 = cell(1,length(folders));
spk_shapes = cell(1,length(folders));
baselineCurrent = cell(1,length(folders));
pertAmp = cell(1,length(folders));
pertArea = cell(1,length(folders));
findV = @(x)x(~cellfun(@isempty,strfind({x.name},'RealNeuron'))).data;
findI = @(x)(x(~cellfun(@isempty,strfind({x.name},'PID'))).data + ...
    x(~cellfun(@isempty,strfind({x.name},'ConstantFromFile'))).data);
findP = @(x)x(~cellfun(@isempty,strfind({x.name},'Waveform'))).data;
for ii = (1:length(folders))
    fprintf(1,'Analizing folder "%s" [%d of %d].\n',folders{ii},ii,length(folders));
    trials_  = {};
    t0_ = [];
    spk_shapes_ = [];
    baselineCurrent_ = [];
    pertAmp_ = [];
    pertArea_ = [];
    files = listFiles(folders{ii},'*.h5','all',{'trash'},{'*_kernel.dat'});
    for file = files'
        [ent,info] = loadH5Trace(file{1});
        [tmptrials,tmpt0,tmpspk_shapes,tmpbaselineCurrent,tmppertAmp,tmppertArea] = splitPRCtracesInTrials( ...
            0:info.dt:info.tend, findV(ent),findI(ent),...
            findP(ent),-35,10,3,3);
        trials_ = [trials_,tmptrials];
        t0_ = [t0_,tmpt0'];
        spk_shapes_ = [spk_shapes_,tmpspk_shapes'];
        baselineCurrent_ = [baselineCurrent_,tmpbaselineCurrent'];
        pertAmp_ = [pertAmp_,tmppertAmp'];
        pertArea_ = [pertArea_,tmppertArea'];
    end
    if REMOVE_OUTLIER_TRIALS
        trial_sd = mean(cellfun(@(x)std(diff(x)),trials_));
        trial_mean = mean(cellfun(@(x)mean(diff(x)),trials_));
        trial_fr = cellfun(@(x)mean(diff(x)),trials_);
        idx_to_remove = find(trial_fr>(trial_mean+trial_sd) | trial_fr<(trial_mean-trial_sd));
        trials_(idx_to_remove) = [];
        t0_(idx_to_remove) = [];
        spk_shapes_(:,idx_to_remove) = [];
        baselineCurrent_(idx_to_remove) = [];
        pertAmp_(idx_to_remove) = [];
        pertArea_(idx_to_remove) = [];
    end
    trials{ii} = trials_;
    t0{ii} = t0_;
    spk_shapes{ii} = spk_shapes_;
    baselineCurrent{ii} = baselineCurrent_;
    pertAmp{ii} = pertAmp_;
    pertArea{ii} = pertArea_;
end
dt = info.dt;

expName = regexp(pwd,'\d{8}[A-Z]\d{2}','match');expName = expName{1};
save('prc_spiketimes.mat','trials','t0','spk_shapes',...
    'baselineCurrent','pertAmp','pertArea','folders',...
    'dt','expName')

%% Plot trials sorted by phase and spike heights and shapes
clear
close all
setFigureDefaults;
colors = [26,91,170;...
    ]./255;
load('prc_spiketimes.mat')
fig_directory = '../figures/rasters';
if ~exist(fig_directory)
    mkdir(fig_directory)
end
for ii = 1:length(trials)
    figure_name = sprintf('raster_%s_%s',expName,folders{ii});
    [phi, dphi, expected_isi, idx]=extractPRC(trials{ii},'unperturbed_mean');
    [~,sort_idx] = sort(phi(:,idx));
    perturbed_isi = find(trials{ii}{1}>0,1)-1;
    predicted_spks = cellfun(@(x)x(find(x>0,1)-1),trials{ii}(sort_idx));
    fig = figure(1);clf;
    set(fig,'color','w')
    % Plot the Raster and predicted isis
    ax(1)=axes();
    plotRastergram(trials{ii}(sort_idx));
    hold on
    for k = -perturbed_isi+1:(length(trials{ii}{1})-perturbed_isi)
        plot(mean(expected_isi)*k+predicted_spks,1:length(predicted_spks),...
            'linewidth',1,'color',colors(1,:))
    end
    axis tight
    xlabel('Spiketimes (ms)')
    % Plot the spike shapes and the mean spike shape.
    ax(2) = axes();
    M = size(spk_shapes{ii},1);
    plot(linspace(0,(M-1)*dt*1e3,M),spk_shapes{ii},...
        'color',[.5,.5,.5],'linewidth',0.5)
    hold on
    plot(linspace(0,(M-1)*dt*1e3,M),mean(spk_shapes{ii}'),...
        'linewidth',1.2,'color',colors(1,:))
    axis tight
    plot([0,0]+max(xlim),[-10,-25],'color','k','linewidth',1)
    plot([max(xlim)-2,max(xlim)],[-25,-25],'color','k','linewidth',1)
    % Plot the predicted versus perturbed isi
    ax(3) = axes();
    perturbed_isis = cellfun(@(x)abs(diff(x(perturbed_isi:perturbed_isi+1))),trials{ii});
    plot(expected_isi,perturbed_isis,'ko','markersize',2.5,'markerfacecolor',[.5,.5,.5])
    hold on
    plot([min(perturbed_isis)-2,max(perturbed_isis)+2],[min(perturbed_isis)-2,max(perturbed_isis)+2],...
        'linewidth',1,'color',colors(1,:))
    axis tight
    grid on
    xlabel('expected interspike interval')
    ylabel('perturbed interspike interval')
    set(ax,'box','off','color','none')
    
    set(ax(3),'position',[0.6,0.1,0.35,0.35],'color','w','ylim',[min(perturbed_isis)-2,max(perturbed_isis)+2],...
        'xlim',[min(perturbed_isis)-2,max(perturbed_isis)+2])
    set(ax(2),'position',[0.6,0.5,0.35,0.45],'ytick',[],'ycolor','w',...
        'xtick',[],'xcolor','w','color','w')
    set(ax(1),'position',[.05,0.1,0.45,0.9],'ytick',[],'ycolor','w')
    caption = sprintf(['Experiment %s. Left - Rastergram sorted by phase.',...
        ' In blue are the expected timestamps according to the overall average ',...
        'firing rate (%3.2fms - %3.1fHz). Top right - Spike shapes. Each gray line is the ',...
        'average spike shape accross that trial. Blue is the overall average. ',...
        'Scale bars are 2ms and 15mV. Bottom right - Perturbed versus expected ',...
        'interspike interval.'],expName,mean(expected_isi),1000./mean(expected_isi));
    set(fig,'paperunits','centimeters','papersize',[15, 15],'paperposition',[0, 0, 15, 15])
    print(fig,'-dpdf',sprintf('%s/%s.pdf',fig_directory,figure_name))
    printFigWithCaption(sprintf('%s/%s.pdf',fig_directory,figure_name),caption)
    movefile(sprintf('%s/%sCaption.pdf',fig_directory,figure_name),sprintf('%s/%s.pdf',fig_directory,figure_name));
end
%%
% PRC summary and peak to baseline
clear all
close all
fig_directory = './';

if ~exist(fig_directory)
    mkdir(fig_directory)
end
setFigureDefaults;
colors = [26,91,170;...
    205,10,32;...
    69,116,40;...
    255,165,0;...
    ]./255;

load('prc_spiketimes.mat')
USE_CORRECTED_METHOD = 0;
prc_x = cell(1,length(trials));
prc_y = cell(1,length(trials));
mean_isi = cell(1,length(trials));
figure_name = sprintf('PRC_summary_%s',expName(expName~=' '));
for ii = 1:length(trials)
    [phi, dphi, expected_isi, idx]=extractPRC(trials{ii},'unperturbed_mean');
    prc_edges = linspace(0,1,100);
    phi = phi./repmat(expected_isi,size(phi,2),1)';
    dphi = dphi./repmat(expected_isi,size(phi,2),1)';
    mean_isi{ii} = mean(expected_isi);
    x = phi(phi>=0 & phi<=1);
    y = dphi(phi>=0 & phi<=1);
    [prc_x{ii},prc_y{ii}] = kernelRegression(x, y, prc_edges);
end

counter=0;
[~,sortidx] = sort(1000./cell2mat(mean_isi));
fig = figure(1);clf
step = 0.07;
ax(1) = axes();
append_to_caption = ['(Firing frequency (Hz), ',...
    'Mean charge delivered to the cell (fC),Number of trials) is'];
for ii = sortidx
    plot(prc_x{ii},prc_y{ii}+counter*step,'color',colors(1,:))
    hold on
    plot(prc_x{ii},zeros(size(prc_x{ii}))+counter*step,'--','color',[.5,.5,.5])
    text(prc_x{ii}(end),prc_y{ii}(end)+counter*step,...
        sprintf('%3.1fHz',1000./mean_isi{ii}),'color',colors(2,:))
    counter=counter+1;
    append_to_caption = sprintf('%s; (%3.1f, %3.2f,%d)',...
        append_to_caption,1000./mean_isi{ii},...
        mean(pertArea{ii})*1.0e3,...
        length(trials{ii}));
end
plot([0,0],[0,0.1],'k')
plot([0,0.1],[0,0],'k')
set(ax(1),'box','off','xcolor','w','ycolor','w','layer','bottom','position',[0.05,0.05,0.85,0.6])
ax(2) = axes();
peak_diff = cellfun(@(x,y)abs(max(y(x<0.5))-max(y(x>0.5))),prc_x(sortidx),prc_y(sortidx));
peak_to_baseline = cellfun(@(x,y)abs(max(y(x<0.5))-max(y(x>0.5)))/(abs(max(y(x<0.5)))+abs(max(y(x>0.5)))),...
    prc_x(sortidx),prc_y(sortidx));
plot(1000./cell2mat(mean_isi(sortidx)),peak_diff,'k--o','markersize',3,'markerfacecolor',colors(1,:))
ylabel('Difference of peaks')
ylim([0,max(ylim)]);
ax(3) = axes();
plot(1000./cell2mat(mean_isi(sortidx)),peak_to_baseline,'ko-.','markersize',3,'markerfacecolor',colors(2,:))
ylabel('Peak to baseline')
xlabel('Firing frequency (Hz)')
set(ax(2),'xaxislocation','top','xcolor','w','layer','bottom','xtick',[])
set(ax(2:3),'box','off','position',[0.1,0.7,0.8,0.25])
set(ax(3),'yaxislocation','right','ycolor',colors(2,:),'ylim',[0,1],'layer','top')
caption = sprintf(['Experiment %s. Firing frequency dependency of the ',...
    'phase response curves. Bottom - Phase advance versus the phase ',...
    'of the perturbation normalized to the mean firing frequency. Each',...
    ' curve denotes the cell at a particular firing frequency. The scale is 0.1  (phase). Top - ',...
    'Peak to baseline and difference of peaks versus the firing frequency of the cell. %s.'],expName,append_to_caption);
set(fig,'paperunits','centimeters','papersize',[15, 20],'paperposition',[0, 0, 15,20])
print(fig,'-dpdf',sprintf('%s/%s.pdf',fig_directory,figure_name))
printFigWithCaption(sprintf('%s/%s.pdf',fig_directory,figure_name),caption)
movefile(sprintf('%s/%sCaption.pdf',fig_directory,figure_name),sprintf('%s/%s.pdf',fig_directory,figure_name));

% PRC and ISI histograms
clear
close all
setFigureDefaults;
colors = [26,91,170;...
    205,10,32;...
    69,116,40;...
    255,165,0;...
    ]./255;
load('prc_spiketimes.mat')
fig_directory = './spk_advance_prc';
if ~exist(fig_directory)
    mkdir(fig_directory)
end
for ii = 1:length(trials)
    
    figure_name = sprintf('spike_advance_prc_%s_%s',expName(expName~=' '),folders{ii});
    [phi, dphi, expected_isi, idx]=extractPRC(trials{ii},'unperturbed_mean');
    [~,sort_idx] = sort(phi(:,idx));
    %perturbed_isi = find(trials{ii}{1}>0,1)-1;
    %predicted_spks = cellfun(@(x)x(find(x>0,1)-1),trials{ii}(sort_idx));
    fig = figure(1);clf;
    set(fig,'color','w')
    % Plot the the histogram of phi
    ax1(1) = axes();
    edges = (min(phi(:))-1:1:max(phi(:))+1);
    counts = histc(phi(:),edges);
    bar(edges,counts,'k')
    axis tight
    % Plot the phase response curve (all phases)
    ax0(1)=axes();
    plot(phi,dphi,'ko','markersize',2,'markerfacecolor',[.5,.5,.5],'markeredgecolor',[.5,.5,.5])
    hold on
    axis tight
    xlabel('Perturbation time (ms)')
    ylabel('Spike advance (ms)')
    set(ax0(1),'position',[0.07,0.1,0.2,0.8])
    set(ax1(1),'position',[0.07,0.75,0.2,0.15],'yaxislocation','right',...
        'ydir','reverse','xtick',[],'xcolor','w')
    set([ax0(1),ax1(1)],'box','off')
    linkaxes([ax0(1),ax1(1)],'x')
    isi = mean(expected_isi);
    x = phi(phi>=0 & phi<=isi);
    y = dphi(phi>=0 & phi<=isi);
    edges = linspace(min(x(:))-0.1,max(x(:))+0.1,20);
    [nX, muY , sY, nSamples] = binSamples(x, y,edges);
    prc_edges = linspace(min(x),max(x),100);
    [prc_x,prc_y] = kernelRegression(x, y, prc_edges);
    % Plot the the histogram of phi
    ax1(2) = axes();
    counts = histc(x(:),edges);
    bar(edges,counts,'k')
    axis tight
    % Plot only phases between zero and one.
    ax0(2)=axes();
    plot(x,y,'ko','markersize',2,'markerfacecolor',[.5,.5,.5],'markeredgecolor',[.5,.5,.5])
    hold on,axis tight
    tmp_phi = phi(:,1:idx-1);
    tmp_dphi = dphi(:,1:idx-1);
    x = tmp_phi(tmp_phi>=0 & tmp_phi<=isi);
    y = tmp_dphi(tmp_phi>=0 & tmp_phi<=isi);
    counts_right = histc(x(:),edges);
    plot(x,y,'ko','markersize',2,'markerfacecolor',colors(1,:),'markeredgecolor',colors(1,:))
    tmp_phi = phi(:,idx+1:length(phi(1,:)));
    tmp_dphi = dphi(:,idx+1:length(phi(1,:)));
    x = tmp_phi(tmp_phi>=0 & tmp_phi<=isi);
    y = tmp_dphi(tmp_phi>=0 & tmp_phi<=isi);
    counts_left = histc(x(:),edges);
    plot(x,y,'ko','markersize',2,'markerfacecolor',colors(2,:),'markeredgecolor',colors(2,:))
    errorbar(nX,muY,sY,'color','k')
    plot(nX,muY,'ko','markerfacecolor',colors(3,:))
    plot(prc_x,prc_y,'color',colors(4,:),'linewidth',1.5)
    xlabel('Perturbation time (ms)')
    ylabel('Spike advance (ms)')
    set(ax0(2),'position',[0.35,0.1,0.3,0.8])
    set(ax1(2),'position',[0.35,0.75,0.3,0.15],'yaxislocation','right',...
        'ydir','reverse','xtick',[],'xcolor','w','xaxislocation','top','layer','bottom')
    linkaxes([ax0(2),ax1(2)],'x')
    set([ax0(2),ax1(2)],'box','off')
    axes(ax1(2));
    hold on
    b = bar(edges,counts_right,'facecolor',colors(1,:));
    bar(edges,counts_left,'facecolor',colors(2,:))
    ax0(3) = axes();
    edges = linspace(min(expected_isi)-0.5,max(expected_isi)+0.5,20);
    counts = histc(expected_isi,edges);
    bar(edges,counts,'k')
    hold on
    plot([mean(expected_isi),mean(expected_isi)],ylim(),...
        'linewidth',1,'color',colors(1,:))
    axis tight
    grid on
    ylabel('counts')
    xlabel('Interspike interval (ms)')
    set([ax0,ax1],'box','off','color','none')
    set(ax0(3),'position',[0.7,0.1,0.2,0.8],'yaxislocation','right')
    caption = sprintf(['Experiment %s. Firing frequency is %3.1fHz [%3.3fms] ',...
        'Mean charge delivered to the cell %3.3f fC, and %d trials. Left - ',...
        'Phase advance versus the phase of the perturbation. All orders ',...
        'are shown. Top is the histogram of phase counts. Middle - Phase ',...
        'advance vs phase of the perturbation. The green points are the ',...
        'mean in each bin. Red and blue are the points drawn from the prcs with order ',...
        '-1 and +1 respectivelly.'],expName,1000./mean(expected_isi),mean(expected_isi),...
    mean(pertArea{ii})*1.0e3,length(trials{ii}));
    set(fig,'paperunits','centimeters','papersize',[20, 10],'paperposition',[0, 0, 20, 10])
    print(fig,'-dpdf',sprintf('%s/%s.pdf',fig_directory,figure_name))
    printFigWithCaption(sprintf('%s/%s.pdf',fig_directory,figure_name),caption)
    movefile(sprintf('%s/%sCaption.pdf',fig_directory,figure_name),sprintf('%s/%s.pdf',fig_directory,figure_name));
end
%%