clc;
clear;
close all;


addpath('./saliency_map/');
addpath('./object_proposal/');
addpath('./test_images/');
addpath('./SpatialPyramid/');
VLFEATROOT = './vlfeat/vlfeat-0.9.21/';  
run([VLFEATROOT 'toolbox/vl_setup']);  % Installing VL_FEAT Toolbox.


imgRoot = './test_images/';
imgRootTmp = './test_images_tmp/';
data_dir = 'data';
resultImgRoot = './result_images/';
dest_sal = './saliency_map/test/'; % directory to contain test images for saliency evaluation.
dest_op = './object_proposal/test_images/'; % directory to contain test images for object proposal evaluation. 

imnames=dir([imgRoot '*' 'jpg']);

if 7==exist(dest_sal,'dir')
    rmdir(dest_sal, 's'); % removing the test folder
end
if 7==exist(dest_op,'dir')
    rmdir(dest_op, 's'); % removing the test folder 
end
if 7==exist(resultImgRoot,'dir')
    rmdir(resultImgRoot, 's'); % removing the result folder 
end
if 7==exist(imgRootTmp,'dir')
    rmdir(imgRootTmp, 's'); % removing the result folder 
end


mkdir(dest_sal); % creating the test folder
mkdir(dest_op); % creating the test folder
mkdir(resultImgRoot); % creating the result folder

%% Copying images to the destination folders

for ii=1:length(imnames)
    imname=[imgRoot imnames(ii).name];

    sal_name = [dest_sal imnames(ii).name];
    op_name = [dest_op imnames(ii).name];
    
    command = ['cp ' imname ' ' sal_name] ;
    system(command);
    
    command = ['cp ' imname ' ' op_name] ;
    system(command);
end

%% ------------------------------ Saliency Map ------------------------------
saliency;


%% --------------------- Object proposal & Localisation -----------------------------

addpath('./object_proposal/matlab');
addpath('./object_proposal/cmex');

%-------------- Global variables (needed for callbacks) --------------
global h imgId files imgDir configFile proposals dense_sift;

%-------------- Input --------------
imgDir = 'test_images';
imgId = 1;
configFile = fullfile('./object_proposal/config', 'rp.mat');
%'config/rp_4segs.mat' to sample from 4 segmentations (slower but higher recall)
%'config/rp.mat' to sample from 1 segmentations (faster but lower recall)

%-------------- Find images in dir: --------------
% files = dir(imgDir);
files = dir([imgRoot '*' '.jpg']);
assert(numel(files) >= 1);
% files = files(3 : end);
if(strcmp(files(1).name, '.svn'))
  files = files(2 : end);
end


%-------------- Processing: --------------
detected = zeros(length(files),1);  % True positives 
% false = zeros(length(files),1);  % False positive

for i = 1:length(files)
    time = tic;
    fprintf('Evaluating Image %d\n', i);
    imgId = i;
    InteractiveCenterDemo(configFile);  % Object Proposals
    imname = files(i).name;
    
    load([imgRoot imname(1:end-4) '.mat'], 'X'); % Load the ground truth bounding box [xmin ymin xmax ymax]
    
    RGB = imread([imgRoot imname]);
    [real_im,w]=removeframe([imgRoot imname]); % run a pre-processing to remove the image frame 
    [m,n,k] = size(real_im);
    command = ['rm ' imgRoot imname(1:end-4) '.bmp'];
    system(command);

    sal_im = imread(['./saliency_map/saliencymap/' imname(1:end-4) '.png']); % saliency map
    spname=['./saliency_map/superpixels/' imname(1:end-4)  '.dat']; 
    superpixels = ReadDAT([m,n],spname); % superpixels matrix
    spnum=max(superpixels(:)); % the actual superpixel number
    
    saliency_contrast = zeros(length(proposals),1);
    area_all = zeros(length(proposals),1);
   
    for r=1:length(proposals)
        box = proposals(r,:); % [ymin, xmin, ymax, xmax]
        area_r = (box(3) - box(1)) * (box(4) - box(2));
        area_all(r) = area_r;
    end
    
    max_area = max(area_all);
    lar_pow_10 = -1;
    while max_area~=0
        lar_pow_10 = lar_pow_10 + 1;
        max_area = floor(max_area/10);
    end
    
    sigma = std(area_all)/(10^lar_pow_10);
    hog_values = zeros(length(proposals),1);
    RS_ = zeros(length(proposals),1);
    NS_ = zeros(length(proposals),1);
    
    for r=1:length(proposals)
        
        box = proposals(r,:); % [ymin, xmin, ymax, xmax]
        area_r = (box(3) - box(1)) * (box(4) - box(2));
        temp = sal_im(box(2):box(4),box(1):box(3));
        sum_of_saliency_region = sum(temp(:))/255;
        RS = (1/area_r)*sum_of_saliency_region; % saliency score w.r.t the proposal region
        RS_(r) = RS;
        
        xmin = box(2)-1;
        ymin = box(1)-1;
        xmax = box(4)+1;
        ymax = box(3)+1;
        if xmin==0
            xmin = 1;
        end
        if ymin==0
            ymin = 1;
        end
        if xmax>=size(superpixels,1)
            xmax=size(superpixels,1);
        end
        if ymax>=size(superpixels,2)
            ymax=size(superpixels,2);
        end
        
        temp = superpixels(xmin:xmax,ymin:ymax);
        temp = unique(temp(:));
        sum_of_saliency_adj = 0.0;
        map = containers.Map(unique(superpixels(:)),ismember(unique(superpixels(:)),temp)); % finding the adjecent superpixels
        
        map = cell2mat(map.values); % converting hashmap to array for O(1) accesses
        area_adj = 0;
        for ii=1:size(superpixels,1)
            for jj=1:size(superpixels,2)
               if map(superpixels(ii,jj))==1
                   if ~((ii>=box(2) && ii<=box(4)) && (jj>=box(1) && jj<=box(3)))
                       sum_of_saliency_adj = sum_of_saliency_adj + (sal_im(ii,jj)/255);
                       area_adj = area_adj + 1;
                   end
               end
            end
        end
        
        NS = (1/area_adj)*sum_of_saliency_adj; % saliency score w.r.t the regions adjecent to the proposal
        NS_(r) = NS;
        
        SC = exp(area_r/sigma.^2)*(RS-NS); % Saliency contrast.
        
        saliency_contrast(r) = SC;
              
    end
   
    idx = find(saliency_contrast==max(saliency_contrast));
%     box_n = proposals(find(area_all==max(area_all(idx))),:); % Selects the largest proposal box. (Not the right way.. Just for the pleasure!!!!!)
    box_n = proposals(idx,:);
    
    spatial_pyramid = zeros(size(idx,1),1000);
    
    if 7==exist(imgRootTmp,'dir')
        rmdir(imgRootTmp, 's'); % removing the result folder 
    end
    mkdir(imgRootTmp);
	
    cccc = 0;
%     
%     %% ---------------temp----------------
%       figure;
%       hold on;
%       imshow(real_im);
%       for i = 1:size(box_n,1)
%         rectangle('Position',[box_n(i,1),box_n(i,2),box_n(i,3)-box_n(i,1),box_n(i,4)-box_n(i,2)],...
%                   'EdgeColor', 'b',...
%                  'LineWidth',1,'LineStyle','-');
%       end
%       axis on;
%       hold off;
%     %% ----------------------------------
    
    for i=1:size(box_n,1)
        outpath = [imgRootTmp '_' num2str(i) '.jpg'];
        if box_n(i,4)-box_n(i,2)<5 && box_n(i,3)-box_n(i,1)>=5
            cccc = cccc+1;
            if box_n(i,2)-5<=0 
                cccc = cccc - 1;
                imwrite(RGB(1:box_n(i,4)+5,box_n(i,1):box_n(i,3),:),outpath);
%             elseif box_n(i,1)-5<=0 && box_n(i,2)-5>0
%                 cccc = cccc - 1;
%                 imwrite(RGB(box_n(i,2)-5:box_n(i,4),1:box_n(i,3)+5,:),outpath);
%             elseif box_n(i,1)-5<=0 && box_n(i,2)-5<=0
%                 cccc = cccc - 1;
%                 imwrite(RGB(1:box_n(i,4)+5,1:box_n(i,3)+5,:),outpath);
            elseif box_n(i,2)-5>0 
                cccc = cccc - 1;
                imwrite(RGB(box_n(i,2)-5:box_n(i,4),box_n(i,1):box_n(i,3),:),outpath);
            
            end     
            
        elseif box_n(i,3)-box_n(i,1)<5 && box_n(i,4)-box_n(i,2)>=5 
            cccc = cccc+1;
%             if box_n(i,2)-5<=0 && box_n(i,1)-5>0
%                 cccc = cccc - 1;
%                 imwrite(RGB(1:box_n(i,4)+5,box_n(i,1):box_n(i,3),:),outpath);
            if box_n(i,1)-5<=0
                cccc = cccc - 1;
                imwrite(RGB(box_n(i,2):box_n(i,4),1:box_n(i,3)+5,:),outpath);
%             elseif box_n(i,1)-5<=0 && box_n(i,2)-5<=0
%                 cccc = cccc - 1;
%                 imwrite(RGB(1:box_n(i,4)+5,1:box_n(i,3)+5,:),outpath);
            elseif box_n(i,1)-5>0
                cccc = cccc - 1;
                imwrite(RGB(box_n(i,2):box_n(i,4),box_n(i,1)-5:box_n(i,3),:),outpath);
            
            end 
            
        elseif box_n(i,4)-box_n(i,2)<5 && box_n(i,3)-box_n(i,1)<5
            cccc = cccc+1;
            if box_n(i,2)-5<=0 && box_n(i,1)-5>0
                cccc = cccc - 1;
                imwrite(RGB(1:box_n(i,4)+5,box_n(i,1):box_n(i,3),:),outpath);
            elseif box_n(i,1)-5<=0 && box_n(i,2)-5>0
                cccc = cccc - 1;
                imwrite(RGB(box_n(i,2)-5:box_n(i,4),1:box_n(i,3)+5,:),outpath);
            elseif box_n(i,1)-5<=0 && box_n(i,2)-5<=0
                cccc = cccc - 1;
                imwrite(RGB(1:box_n(i,4)+5,1:box_n(i,3)+5,:),outpath);
            elseif box_n(i,1)-5>0 && box_n(i,2)-5>0
                cccc = cccc - 1;
                imwrite(RGB(box_n(i,2)-5:box_n(i,4),box_n(i,1)-5:box_n(i,3),:),outpath);
            
            end        
                
        elseif box_n(i,4)-box_n(i,2)>=5 && box_n(i,3)-box_n(i,1)>=5         
            imwrite(RGB(box_n(i,2):box_n(i,4),box_n(i,1):box_n(i,3),:),outpath);
        end
    end
    
    imnames=dir([imgRootTmp '*' '.jpg']);
    
    for i=1:size(box_n,1)
        
%         % ------ Calculation of Dense SIFT desriptor ------
%         [f1,d1] = vl_dsift(single(rgb2gray(real_im)),'bounds',[box_n(i,2),box_n(i,1),box_n(i,4),box_n(i,3)], 'size', 10, 'step', 20);
%         kd = KDTreeSearcher(double(d1)');
%         s = struct('features',f1,'descr',d1,'KDTree',kd);
%         if i==1
%             dense_sift = s;
%         else            
%             dense_sift = [dense_sift;s];
%         end
        
        % ------- Spatial Pyramid Matching --------
        filename = cell(1,1); 
        filename{1} = imnames(i).name;
        params.maxImageSize = 1000
        params.gridSpacing = 1
        params.patchSize = 16
        params.dictionarySize = 200
        params.numTextonImages = 50
        params.pyramidLevels = 2
        pyramid = BuildPyramid(filename,imgRootTmp,data_dir,params);
        spatial_pyramid(i,:) = pyramid;                       
        
    end
    axis on;
    
    kd = KDTreeSearcher(spatial_pyramid); % Converting to K-d tree, inorder to speed up the KNN search.
    tmpvar = 11;
    if size(kd.X,1)<11
        tmpvar = size(kd.X,1);
    end
    groups = zeros(size(spatial_pyramid,1),tmpvar);
    scores = zeros(size(spatial_pyramid,1),1);
    for i=1:size(spatial_pyramid,1)
        groups(i,:) = knnsearch(kd,spatial_pyramid(i,:),'K',tmpvar);    
        score = 0.0;
        for j=1:tmpvar
            for k=j+1:tmpvar
               score = score + norm(spatial_pyramid(groups(i,j),:)' - spatial_pyramid(groups(i,k),:)');
            end
        end
        scores(i) = score;
    end
    
    scores_sorted = sort(scores,'ascend');
    tmpvar = 5;
    if length(scores_sorted)<tmpvar
        tmpvar = length(scores_sorted);
    end
    
    scores_sorted = scores_sorted(1:tmpvar);
    
    
%     x_min = 0;
%     x_max = 0;
%     y_min = 0;
%     y_max = 0;
    x_min = Inf;
    y_min = Inf;
    y_max = 0;
    x_max = 0;
    count = 0;
    for i=1:length(scores_sorted)
        
        idx_ = find(scores==scores_sorted(i));
        count = count + length(idx_);
        q = mod(i,8);
        if q==0
            q = 1;
        end
%         sample = [box_n(idx_,1),box_n(idx_,2),box_n(idx_,3)-box_n(idx_,1),box_n(idx_,4)-box_n(idx_,2)];
%         rectangle('Position',sample(1,:),...                
%                    'EdgeColor', edge_colors(q),...  
%                    'LineWidth',1,'LineStyle','-');
%         y_min = y_min + sum(box_n(idx_,1));
%         y_max = y_max + sum(box_n(idx_,3));
%         x_min = x_min + sum(box_n(idx_,2));
%         x_max = x_max + sum(box_n(idx_,4));

        y_min = min([y_min box_n(idx_,1)']);
        x_min = min([x_min box_n(idx_,2)']);
        y_max = max([y_max box_n(idx_,3)']);
        x_max = max([x_max box_n(idx_,4)']);
      
    end
    
%     y_min = ceil(y_min/count);
%     y_max = ceil(y_max/count);
%     x_min = ceil(x_min/count);
%     x_max = ceil(x_max/count);
    
    if isinf(y_min) || isnan(y_min) || y_min == 0
        y_min = 1;
    end    
    if isinf(x_min) || isnan(x_min) || x_min == 0
        x_min = 1;
    end    
    if isinf(y_max) || isnan(y_max) || y_max == 0
        y_max = 2;
    end
    if isinf(x_max) || isnan(x_max) || x_max == 0
        x_max = 2;
    end
    
    fig = figure;
    hold on;
    imshow(real_im);
    edge_colors = ['r','g','b','y','m','c','w','k'];
    rectangle('Position',[y_min,x_min,y_max-y_min,x_max-x_min],...                
               'EdgeColor', edge_colors(2),...  
               'LineWidth',1,'LineStyle','-');
    axis on;
    
    frm = getframe(fig);
    imwrite(frm.cdata, [resultImgRoot imname(1:end-4) '_result.png']); % Storing the result image.
    
    
    %% -------- Evaluating the Localisation ---------
    bbgt = X; % ground truth
    bb = [x_min y_min x_max y_max]; % my answer  
    bi = [max(bb(1),bbgt(1)) ; max(bb(2),bbgt(2)) ; min(bb(3),bbgt(3)) ; min(bb(4),bbgt(4))]; % intersection
    minoverlap = 0.5; % min overlap required for the localisation to be correct ( CorLoc ).
    
    iw=bi(3)-bi(1)+1; % intersection width
    ih=bi(4)-bi(2)+1; % intersection height
    ov=0.0;
    if iw>0 && ih>0                
        % compute overlap as area of intersection / area of union
        % ua -- union area
        % ov -- area of intersection / area of union
        ua=(bb(3)-bb(1)+1)*(bb(4)-bb(2)+1)+...
           (bbgt(3)-bbgt(1)+1)*(bbgt(4)-bbgt(2)+1)-...
           iw*ih;
        ov=iw*ih/ua;
    end
    
    % assign detection as true positive/don't care/false positive
    if ov>=minoverlap
        detected(i) = 1;                 % true positive
%     else
%         false(i)=1;                    % false positive
    end
   
    fprintf('done in %0.2f seconds!\n-------------------------------------\n',toc(time));
    close all;
end








