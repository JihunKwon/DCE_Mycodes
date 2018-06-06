%Gaussian Filter
clear;
clc;

% 1. Change the directory
% 2. Change base folder
% 3. Change PatientID
% 4. Change PatientName

% G1
 basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group1_GNP+IR/Day6_post/16-BERB035_334_334_20160412__E5_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group1_GNP+IR/Day6_post/16-BERB035_337_337_20160412__E5_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group1_GNP+IR/Day6_post/16-BERB035_340_340_20160412__E5_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group1_GNP+IR/Day6_post/16-BERB035_345_345_20160412__E4_P1';
% G2
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group2_IR/Day6_post/16-BERB035_351_351_20160419_post__E5_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group2_IR/Day6_post/16-BERB035_352_352_20160419_post__E3_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group2_IR/Day6_post/16-BERB035_353_353_20160419_post__E3_P1';
% G3
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group3_GNP/Day6_post/BERB035_338_338_20170503_v2__E3_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group3_GNP/Day6_post/BERB035_343_343_20170503_v2__E3_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group3_GNP/Day6_post/BERB035_344_344_20170503_v3__E3_P1';
% G4
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group4_Control/Day6_post/16-BERB035_342_342_20160426_v2__E3_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group4_Control/Day6_post/16-BERB035_350_350_20160426_v2__E3_P1';
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group4_Control/Day6_post/16-BERB035_355_355_20160426_v2__E3_P1';

%% Apply for G-filtered images
% basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Group1_GNP+IR/Day6_post/New_337/T1G1_337';

for i=1:450
    Y=0;
    fdcm_pre = sprintf('MRIm%03d.dcm', i);
    fname = fullfile(basefolder, fdcm_pre);
    [X,map] = dicomread(fname);
    info = dicominfo(fname);
    
%     %Get next frame
%     if i <= 445
%         fdcm_pre2 = sprintf('MRIm%03d.dcm', i+5);
%     else
%         fdcm_pre2 = sprintf('MRIm%03d.dcm', i); % last frame is not averaged
%     end
%     fname2 = fullfile(basefolder, fdcm_pre2);
%     [X2,map2] = dicomread(fname2);
%     info2 = dicominfo(fname2);
%     
%     %Get next next frame
%     if i <= 440
%         fdcm_pre3 = sprintf('MRIm%03d.dcm', i+10);
%     else
%         fdcm_pre3 = sprintf('MRIm%03d.dcm', i);
%     end
%     fname3 = fullfile(basefolder, fdcm_pre3);
%     [X3,map3] = dicomread(fname3);
%     info3 = dicominfo(fname3);
%     
%     %% For T2G0 (Working well. Useful filter)
%     info.PatientID = '334_T2G0';
%     info.PatientName.FamilyName = '334_T2G0';
%     Y2_0 = (X+X2)/2;
%     fdcm_2 = sprintf('2_0_MRIm%03d.dcm', i);
%     dicomwrite(Y2_0, fdcm_2, info);
%     
%      %% For T3G0
%     info.PatientID = '334_T3G0';
%     info.PatientName.FamilyName = '334_T3G0';
%     Y3 = (X+X2+X3)/3;
%     fdcm_3 = sprintf('3_0_MRIm%03d.dcm', i);
%     dicomwrite(Y3, fdcm_3, info);

    %Get next frame
    if i > 5
        fdcm_pre2 = sprintf('MRIm%03d.dcm', i-5);
    else
        fdcm_pre2 = sprintf('MRIm%03d.dcm', i); % last frame is not averaged
    end
    fname2 = fullfile(basefolder, fdcm_pre2);
    [X2,map2] = dicomread(fname2);
    info2 = dicominfo(fname2);
    
    %% For T2G0 (Working well. Useful filter)
    info.PatientID = '334_T2G0';
    info.PatientName.FamilyName = '334_T2G0';
    Y2_0 = (X+X2)/2;
    fdcm_2 = sprintf('2_0_MRIm%03d.dcm', i);
    dicomwrite(Y2_0, fdcm_2, info);
end