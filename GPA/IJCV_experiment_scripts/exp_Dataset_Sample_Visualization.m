clear all
close all
clc


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'}
samplePick = [1, 1, 1, 1, 1, 1];
markerSize = [6, 6, 6, 6, 6, 3] ;

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

addpath('../DefGPA', '../');


figSize = [100, 50, 250, 200];

figformat = 'pdf';

rescaleFactor = [0.4, 0.25, 0.25, 0.4, 0.4, 0.4];

for dcnt = [1, 2, 3]
    
      datasetName = datasetCell{dcnt};

     [DataCell, VisCell, SampleCell] = readDataset (datasetName);
     
     sampleCnt = samplePick(dcnt);
     
     markerData = rescaleFactor(dcnt) * DataCell{sampleCnt};
     
     sampleName = SampleCell{sampleCnt};
     
     hFig = figure(dcnt);
     
     img = imread(sampleName);
     if(dcnt == 1)
         markerData = rescaleFactor(dcnt) * (DataCell{sampleCnt} - [1; 20]);
         img = imcrop(img, [1, 20, 720-1, 540-1]);
     end
     
     img = imresize(img, rescaleFactor(dcnt));
     
     imshow(img);
     
     hold on;
     
     plot(markerData(1,:), markerData(2,:), 'g.', 'MarkerSize', markerSize(dcnt));
     
     hold off;
     
     axis off
     
     ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
     
     set(hFig, 'Position', figSize);
     
     fig_saved_name = [fdir,  'exp_dataset_sample_', datasetName];
     
     ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','vector');
     
     pause(0.5);
     
end



for dcnt = 4
        
      datasetName = datasetCell{dcnt};

     [DataCell, VisCell, Samples] = readDataset (datasetName);
     
     sample1 = Samples{1};
     sample2 = Samples{2};
     sample3 = Samples{3};
     sample4 = Samples{4};
          
     img1 = imread(sample1);
     img2 = imread(sample4);
     
    hFig = figure(dcnt);

    subplot('Position', [0.04, 0.05, 0.45, 0.9]);
    imshow(img1);
    
    subplot('Position', [0.04+0.0+0.45, 0.05, 0.45, 0.9]);
    imshow(img2);
     
     ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
     
     set(hFig, 'Position', figSize);
     
     fig_saved_name = [fdir,  'exp_dataset_sample_', datasetName];
     
     exportgraphics(hFig, [fig_saved_name, '.', figformat], 'Resolution', 300, 'BackgroundColor','k');
     
     pause(0.5);
     
    
end




for dcnt = 5
    
      datasetName = datasetCell{dcnt};

     [DataCell, VisCell, SampleCell] = readDataset (datasetName);
     
     sampleCnt = samplePick(dcnt);
     
     markerData = DataCell{sampleCnt};
     
     sampleName = SampleCell{sampleCnt};
     
     hFig = figure(dcnt);
     
     img = imread(sampleName);
     
     img = imresize(img, rescaleFactor(dcnt));
     
     imshow(img);

     axis off
     
     ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
     
     set(hFig, 'Position', figSize);
     
     fig_saved_name = [fdir,  'exp_dataset_sample_', datasetName];
     
     ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','vector');
     
     pause(0.5);
    
end



for dcnt = 6
    
      datasetName = datasetCell{dcnt};

     [DataCell, VisCell, SampleCell] = readDataset (datasetName);
     
     sampleCnt = samplePick(dcnt);

     [size1, size2, ~] = size(imread(SampleCell{sampleCnt}{1}));
   
     hFig = visualizeGenenerativelModel (dcnt, sampleCnt, markerSize(dcnt), rescaleFactor(dcnt));
     
     figSize = [100, 50, 500, 220];
     
     set(hFig, 'Position', figSize);
        
     fig_saved_name = [fdir,  'exp_dataset_sample_', datasetName];
     
     exportgraphics(hFig, [fig_saved_name, '.', figformat], 'ContentType','vector');
     
     pause(0.5);
    
end










function hFig = visualizeGenenerativelModel (fid, sampleChoice, markersize, rescaleFactor)

allSamples = 1 : 3 : 600;

fdir = '../../dataset/stereoToyRug/';

% number of frames
n = 651;

% number of points
m = 30;

% read left image point tracks
Xl = load([fdir, 'X_left.txt']);

% read right image point tracks
Xr = load([fdir, 'X_right.txt']);

% read 3D points
X = load([fdir, 'X.txt']);

% read camera parameters and reassemble them to projection matrix
Kl = load([fdir, 'Kl.txt']);
Kr = load([fdir, 'Kr.txt']);
Rl = load([fdir, 'Rl.txt']);
Rr = load([fdir, 'Rr.txt']);
Tl = load([fdir, 'Tl.txt']);
Tr = load([fdir, 'Tr.txt']);
Pl = Kl*[Rl  100*Tl];
Pr = Kr*[Rr 100*Tr];


% visualise data sequentially
hFig = figure(fid);


for i = allSamples (sampleChoice)
    
    fprintf('frame %03d / %03d\n',i,n);
    
    % read left image
    Il = imread(sprintf([fdir, 'I_left_%03d.jpg'],i));
    Il = imresize(Il, rescaleFactor);
    
    % read right image
    Ir = imread(sprintf([fdir, 'I_right_%03d.jpg'],i));
    Ir = imresize(Ir, rescaleFactor);
    
    % extract left points
    ql = rescaleFactor * Xl(2*i-1:2*i,:);

    % extract right points
    qr = rescaleFactor * Xr(2*i-1:2*i,:);

    % extract 3D points
    Q = 100 * X(3*i-2:3*i,:);

    % reproject 3D points
    hql = rescaleFactor * Pl*[Q ; ones(1,m)];
    hql = rescaleFactor * hql(1:2,:)./repmat(hql(3,:),2,1);
    hqr = rescaleFactor * Pr*[Q ; ones(1,m)];
    hqr = rescaleFactor * hqr(1:2,:)./repmat(hqr(3,:),2,1);
    

    subplot('Position', [0.04, 0.05, 0.45, 0.9]);
    hold off, imshow(Il), hold on
    plot(ql(1,:),ql(2,:),'go','markerfacecolor','g', 'MarkerSize', markersize);
    plot(hql(1,:),hql(2,:),'bx',  'MarkerSize', markersize);
    ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1; ax.Title.FontSize = 12;
    title('Left');
    

    subplot('Position', [0.04+0.02+0.45, 0.05, 0.45, 0.9]);
    hold off, imshow(Ir), hold on
    plot(qr(1,:),qr(2,:),'go','markerfacecolor','g', 'MarkerSize', markersize);
    plot(hqr(1,:),hqr(2,:),'bx',  'MarkerSize', markersize);
    ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1; ax.Title.FontSize = 12;
    title('Right');
    

    pause(0.1);
    
end

end





