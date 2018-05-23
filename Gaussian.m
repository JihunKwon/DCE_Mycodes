%Gaussian Filter
clear;
clc;

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

for i=1:5
    fdcm_pre = sprintf('MRIm%03d.dcm', i);
    fname = fullfile(basefolder, fdcm_pre);
    [X,map] = dicomread(fname);
    info = dicominfo(fname);
    
    %% For T1G2
    info.PatientID = '334_T1G2';
    info.PatientName.FamilyName = '334_T1G2';
    Y2 = imgaussfilt(X,2);
    fdcm_2 = sprintf('2_MRIm%03d.dcm', i);
    fdcm_2 = char(fdcm_2);
    dicomwrite(Y2, fdcm_2, info);

    %% For T1G4
    info.PatientID = '334_T1G4';
    info.PatientName.FamilyName = '334_T1G4';
    Y4 = imgaussfilt(X,4);
    fdcm_4 = sprintf('4_MRIm%03d.dcm', i);
    fdcm_4 = char(fdcm_4);
    dicomwrite(Y2, fdcm_4, info);
end