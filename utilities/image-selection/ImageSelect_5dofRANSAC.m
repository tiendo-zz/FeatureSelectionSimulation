function [PairWiseMatches,PairWiseMatches_homo,featureExtracted,match_id,ImageId,t] = ImageSelect_5dofRANSAC(img_idx,featureExtracted,match_id_candid,intrinsic,measurePosesId,rate,NumOfPoses,t)
% img_idx : the new image's index
% VT_flag : the index of image which should be used to do VT matching
% match_id : if the new image won't be sent to do VT matching, then match
% it with the images inside match_id 
global datapath image_folder CameraModel output_folder config_folder

%% Extratc FREAK features for current image
dlmwrite([output_folder,'/match_id_',num2str(img_idx),'.txt'],sort(measurePosesId(match_id_candid)));
command = [datapath,'ExtractFREAK.sh ',image_folder,' ',config_folder,' ',output_folder,' ',num2str(rate),' ',num2str(img_idx),' ',num2str(0)];
tic;
system(command);
t(img_idx,1) = tic;
ImageId = load([output_folder,'/idx_selected.txt']);
if ImageId > NumOfPoses
    match_id = [];
    PairWiseMatches = cell(1,1);
    PairWiseMatches_homo = cell(1,1);
    featureExtracted= cell(1,1);
    return;
end


new_FREAK = load([output_folder,'/features_clone/',num2str(ImageId),'.txt']);
featureExtracted_new = new_FREAK(2:end,1:2)';
featureExtracted{img_idx} = featureExtracted_new;


% run query VT matching and update vt tree
% return number of matches & match image idx & features
command = [datapath,'QueryVT.sh ',image_folder,' ',config_folder,' ',output_folder,' ',num2str(ImageId),' ',num2str(0)];
system(command);
% check if get matching images
match_id_vt = load([output_folder,'/match_data/MatchId_',num2str(ImageId),'.txt']);
match_id = match_id_candid;
if ~isempty(match_id_vt)
   match_id_vt = match_id_vt(find(match_id_vt(:,2)>10),:);
   for i = 1:size(match_id_vt,1)
       if isempty(find(match_id == match_id_vt(i,1))) % this is a new match image 
          match_id = [match_id_vt(i,1) match_id];
       end
   end
else
   match_id = match_id_candid;
end


PairWiseMatches = cell(nchoosek(img_idx,2),1);
PairWiseMatches_homo = cell(nchoosek(img_idx,2),1);

if ~isempty(match_id)

   for i=1:length(match_id)
       PairWiseMatches{SUTrilInd(img_idx,match_id(i),img_idx)} = load([output_folder,'/match_data/5doefInlier_',num2str(measurePosesId(match_id(i))),'_',num2str(ImageId),'.txt']);
       PairWiseMatches_homo{SUTrilInd(img_idx,match_id(i),img_idx)} =  PairWiseMatches{SUTrilInd(img_idx,match_id(i),img_idx)};
    
       if isempty(PairWiseMatches{SUTrilInd(img_idx,match_id(i),img_idx)})
          continue;
       end
    
       homo_1 = StaticUndistort(intrinsic.fc,intrinsic.cc,intrinsic.kc,PairWiseMatches{SUTrilInd(img_idx,match_id(i),img_idx)}(:,1:2)',CameraModel);
       homo_2 = StaticUndistort(intrinsic.fc,intrinsic.cc,intrinsic.kc,PairWiseMatches{SUTrilInd(img_idx,match_id(i),img_idx)}(:,4:5)',CameraModel);    
       PairWiseMatches_homo{SUTrilInd(img_idx,match_id(i),img_idx)}(:,1:2) = homo_1';
       PairWiseMatches_homo{SUTrilInd(img_idx,match_id(i),img_idx)}(:,4:5) = homo_2';
   end
end
