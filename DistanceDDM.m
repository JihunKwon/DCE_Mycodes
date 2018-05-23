clear;
clc;
%% Prepare figure
figure(11); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(12); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(13); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(14); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(21); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(22); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(23); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(24); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(31); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(32); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(33); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(34); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(41); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(42); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(43); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(44); clf('reset');set(gcf,'Position',[440   378   560   420]);

%%
%Convert voxel to volume
VoxelVolume = 0.15626 * 0.15625 * 1; % [mm]
basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Analysis/Input';

AnimalID = [334 337 338 340 342 343 344 345 350 351 352 353 355];
AnimalID = AnimalID(:);
%1:334, 2:337, 3:338, 4:340, 5:342, 6:343, 7:344, 8:345, 9:350,
%10:351, 11:352, 12:353, 13:355
%GdNP+RT:1,2,4,8; RT:10,11,12; GdNP:3,6,7; Control:5,9,13

DepthBinLimit = 70;
TotCount = zeros(13,5,DepthBinLimit);
Total = zeros(13,5,DepthBinLimit);
NetCount = zeros(13,5,DepthBinLimit);
Net = zeros(13,5,DepthBinLimit);
Ave = zeros(13,5,DepthBinLimit);

NetCountBin = zeros(13,5,20);
NetBin = zeros(13,5,20);
AveBin = zeros(13,5,20);
NetCountBinMax = zeros(13,20);
NetBinMax = zeros(13,20);
AveBinMax = zeros(13,20);

NetCountBinSiju = zeros(13,5,20);
NetBinSiju = zeros(13,5,20);
AveBinSiju = zeros(13,5,20);
NetCountBinMaxSiju = zeros(13,20);
NetBinMaxSiju = zeros(13,20);
AveBinMaxSiju = zeros(13,20);

TargetLog = zeros(13,1);
KtransFrac = zeros(13,DepthBinLimit);
AveKtrans = zeros(13,1);

NetCountBinRelative = zeros(13,5,20);
NetBinRelative = zeros(13,5,20);
AveBinRelative = zeros(13,5,20);
NetCountBinMaxRelative = zeros(13,20);
NetBinMaxRelative = zeros(13,20);
AveBinMaxRelative = zeros(13,20);

MaxDepth = zeros(13,6); %Relative depth from edge. (AnimalID, Slice, MaxDepth)
RelativeCount = zeros(13,5);
RelativeNet = zeros(13,5);
RelativeAve = zeros(13,5);
%% Main Program
for k = 1:size(AnimalID,1)
    %Import Tumor mask
    fname_mask = sprintf('%d - as a 90 frames MultiVolume by ImagePositionPatientInstanceNumber frame 89-label.nrrd', AnimalID(k));
    fname = fullfile(basefolder, fname_mask);
    
    if exist(fname_mask,'file')
        [Vmask,mask] = nrrdread2(fname);
        %imshow(Vmask(:,:,2),[]);
         
        %z=1 is the lowest slice, z=5 is the highest.
%        for z = 1:5 %Comment this if only maximum slice is needed.
%             %**Calculate all slices**
%             TargetSlice = z; 
           %**Or find the largest tumor slice and escape for loop**
            z=1;
%             TargetSlice = maxslice(Vmask); %In case of target==MaxVolume
%             TargetLog(k,1) = TargetSlice;  %In case of target==MaxVolume
            TargetSlice = 3; %In case of target==center (slice #3)
            
            Vmask_target = Vmask(:,:,TargetSlice); %Target slice into 2D
            Vddm = bwdist(~Vmask_target);%Euclidean Distance Transformation
            %imshow(Vddm(:,:),[]);

            Vddm_distance = Vddm*0.15624999999999997;%Convert DDM to real distance [mm]
            %imshow(Vddm_distance(:,:),[]);
            %imcontour(Vddm_distance,7);

            %Import Ktrans map
            fktrans = sprintf('%d_Output Ktrans image_pop.nrrd', AnimalID(k));
            fname_ktrans = fullfile(basefolder, fktrans);

            if exist(fname_ktrans,'file')
                [Vktrans,ktrans] = nrrdread2(fname_ktrans);
                Vktrans_target = Vktrans(:,:,TargetSlice);
                %imshow(Vktrans_target(:,:),[]);
                AveKtrans(k,1) = sum(Vktrans(:,:,3))/sum(Vmask(:,:,3));
                %% Distance-corelated analysis
                MaxDDM = 3.4375; %This has to be carefully defined.
                DDMBIN = 15;
                for i = 1:1:MaxDDM*DDMBIN+1
                    DDM = Vddm_distance*DDMBIN;
                    %Gradually shrink DDM
                    DDM(i > DDM) = 0;
                    DDM(i <= DDM) = 1;
                    %imshow(DDM(:,:),[]);

                    %Count voxel number
                    %Since Ktrans value sometimes becomes zero, voxel number needs to be
                    %calculated before appying DDM to map
                    Value = (DDM ~= 0.0);
                    TotCount(k,z,i) = sum(Value(:)); 

                    %Apply DDM to Ktrans map
                    DDMmask = DDM.*Vktrans_target; 

                    %Calculate sum and voxel number of each mask
                    Total(k,z,i) = sum(DDMmask(:)); %sum
                end

                %To include all the voxels, first bin needs to be one.
                for i = 1:1:MaxDDM*DDMBIN 
                    Net(k,z,i) = Total(k,z,i) - Total(k,z,i+1);
                    NetCount(k,z,i) = TotCount(k,z,i) - TotCount(k,z,i+1);
                    Ave(k,z,i) = Net(k,z,i)/NetCount(k,z,i);
                end
                
                %% Absolute Depth
                %If bin=1, spectrum is too noisy. Modulate spectrum by "ModulateBin".
                a = 1; 
                b = 1;
                ModulateBin = 5;
                ModulateBinSiju = 10;
                for i = 1:1:MaxDDM*DDMBIN
                    if ((a-1)*ModulateBin<i && i<=a*ModulateBin)
                        NetCountBin(k,z,a) = NetCountBin(k,z,a) + NetCount(k,z,i);
                        NetBin(k,z,a) = NetBin(k,z,a) + Net(k,z,i);
                        NetCountBinMax(k,a) = squeeze(NetCountBin(k,z,a));
                        NetBinMax(k,a) = squeeze(NetBin(k,z,a));
                    end
                    if rem(i,ModulateBin) == 0
                        AveBin(k,z,a) = NetBin(k,z,a)/NetCountBin(k,z,a);
                        AveBinMax(k,a) = squeeze(AveBin(k,z,a));
                        a = a+1;
                    end
                end
                for i = 1:1:MaxDDM*DDMBIN
                    if ((b-1)*ModulateBinSiju<i && i<=b*ModulateBinSiju)
                        NetCountBinSiju(k,z,b) = NetCountBinSiju(k,z,b) + NetCount(k,z,i);
                        NetBinSiju(k,z,b) = NetBinSiju(k,z,b) + Net(k,z,i);
                        NetCountBinMaxSiju(k,b) = squeeze(NetCountBinSiju(k,z,b));
                        NetBinMaxSiju(k,b) = squeeze(NetBinSiju(k,z,b));
                    end
                    if rem(i,ModulateBinSiju) == 0
                        AveBinSiju(k,z,b) = NetBinSiju(k,z,b)/NetCountBinSiju(k,z,b);
                        AveBinMaxSiju(k,b) = squeeze(AveBinSiju(k,z,b));
                        b = b+1;
                    end
                end
                
                %% Relative Depth
                for i=DepthBinLimit:-1:1
                    depth = FindMaxDepth(Net,k,i); %Find Maximum i. i=depth
                    if depth~=0
                        MaxDepth(k,1) = depth;
                        break;
                    end
                end
                
                for j = 2:6
                    MaxDepth(k,j) = MaxDepth(k,1)/100*(j-1)*20; %Convert to 20%,40%,...,100%
                end
                MaxDepth(k,1)=0; %To use this variable in the next for loop
                
                for i = 1:1:DepthBinLimit
                    if i <= MaxDepth(k,2)
                        RelativeCount(k,1) = RelativeCount(k,1) + NetCount(k,1,i);
                        RelativeNet(k,1) = RelativeNet(k,1) + Net(k,1,i);
                    elseif i <= MaxDepth(k,3)
                        RelativeCount(k,2) = RelativeCount(k,2) + NetCount(k,1,i);
                        RelativeNet(k,2) = RelativeNet(k,2) + Net(k,1,i);
                    elseif i <= MaxDepth(k,4)
                        RelativeCount(k,3) = RelativeCount(k,3) + NetCount(k,1,i);
                        RelativeNet(k,3) = RelativeNet(k,3) + Net(k,1,i);
                    elseif i <= MaxDepth(k,5)
                        RelativeCount(k,4) = RelativeCount(k,4) + NetCount(k,1,i);
                        RelativeNet(k,4) = RelativeNet(k,4) + Net(k,1,i);
                    elseif i <= MaxDepth(k,6)
                        RelativeCount(k,5) = RelativeCount(k,5) + NetCount(k,1,i);
                        RelativeNet(k,5) = RelativeNet(k,5) + Net(k,1,i);
                    end
                end
                
                for i = 1:1:DepthBinLimit
                    for j = 1:5
                        RelativeAve(k,j) = RelativeNet(k,j)/RelativeCount(k,j);
                    end
                end
            else
                fprintf('Animal_%d, Ktrans map not imported\n', AnimalID(k));
            end
            AveBinMax(isnan(AveBinMax))=0;
            AveBinMaxSiju(isnan(AveBinMaxSiju))=0;
%        end %Comment this if only maximum slice is needed.
    else
        fprintf('Animal_%d, mask not imported\n', AnimalID(k));
    end
end

%% Average Ktrans for entire volume
fprintf('Average Ktrans for GdNP+RT = %d\n', (AveKtrans(1,1)+AveKtrans(2,1)+AveKtrans(4,1)+AveKtrans(8,1))/4);
fprintf('Average Ktrans for RT = %d\n', (AveKtrans(10,1)+AveKtrans(11,1)+AveKtrans(12,1))/3);
fprintf('Average Ktrans for GdNP = %d\n', (AveKtrans(3,1)+AveKtrans(6,1)+AveKtrans(7,1))/3);
fprintf('Average Ktrans for Control = %d\n', (AveKtrans(5,1)+AveKtrans(9,1)+AveKtrans(13,1))/3);
                
%% Draw Figure
MaxBin = 60;
RoundMax = round(MaxBin/ModulateBin);
RoundMaxSiju = round(MaxBin/ModulateBinSiju);
for k = 1:size(AnimalID,1)
    if (k == 1 || k == 2 || k == 4 || k == 8) %GdNP+RT
        if k == 1
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 2
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 4
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 8
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 10 || k == 11 || k == 12) %RT
        if k == 10
            figure(21); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(22); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(23); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 11
            figure(21); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(22); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(23); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 12
            figure(21); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(22); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(23); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 3 || k == 6 || k == 7) %GdNP
        if k == 3
            figure(31); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(32); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(33); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 6
            figure(31); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(32); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(33); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 7
            figure(31); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(32); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(33); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 5 || k == 9 || k == 13) %Control
        if k == 5
            figure(41); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(42); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(43); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-p','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 9
            figure(41); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(42); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(43); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 13
            figure(41); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(42); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(43); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,AveBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end
    else
        fprintf('Animal_%d, mask not imported\n', AnimalID(k));
    end
end


%% Add labels (max slice)
%{
%GdNP+RT starts%
figure(42);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Total K^{trans}, GdNP+RT');
xlim([0 50])
ylim([0 820]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Total K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransSum_GdNP+RT_max.eps', '-depsc2');

figure(43);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, GdNP+RT');
xlim([0 50])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP+RT_max.eps', '-depsc2');

figure(13);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, GdNP+RT');
xlim([0 50])
ylim([0 3.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans}');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_max.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(21);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Total K^{trans}, RT');
xlim([0 50])
ylim([0 820]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Total K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransSum_RT_max.eps', '-depsc2');

figure(22);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, RT');
xlim([0 50])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_RT_max.eps', '-depsc2');

figure(23);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, RT');
xlim([0 50])
ylim([0 3.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans}');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_max.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(31);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Total K^{trans}, GdNP');
xlim([0 50])
ylim([0 820]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Total K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransSum_GdNP_max.eps', '-depsc2');

figure(32);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, GdNP');
xlim([0 50])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP_max.eps', '-depsc2');

figure(33);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, GdNP');
xlim([0 50])
ylim([0 3.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans}');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_max.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(41);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Total K^{trans}, Control');
xlim([0 50])
ylim([0 820]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Total K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransSum_Cont_max.eps', '-depsc2');

figure(42);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, Control');
xlim([0 55]);
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_Cont_max.eps', '-depsc2');

figure(43);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, Control');
xlim([0 50]);
ylim([0 3.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans}');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Cont_max.eps', '-depsc2');
%Control ends%
%}


%% Add labels (center slice) (AverageUnderAIFMask mode) 
%{
%GdNP+RT starts%
figure(11);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin), GdNP+RT');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_GdNP+RT_center.eps', '-depsc2');

figure(12);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, GdNP+RT');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP+RT_center.eps', '-depsc2');

figure(13);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, GdNP+RT');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_center.eps', '-depsc2');

figure(14);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, Relative Depth, GdNP+RT');
xlim([0 5])
ylim([0 1.5]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_center_Relative.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(21);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin), RT');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_RT_center.eps', '-depsc2');

figure(22);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, RT');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_RT_center.eps', '-depsc2');

figure(23);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin), RT');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center.eps', '-depsc2');

figure(24);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, Relative Depth, RT');
xlim([0 5])
ylim([0 1.5]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center_Relative.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(31);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin), GdNP');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_GdNP_center.eps', '-depsc2');

figure(32);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, GdNP');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP_center.eps', '-depsc2');

figure(33);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, GdNP');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center.eps', '-depsc2');

figure(34);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, Relative Depth, GdNP');
xlim([0 5])
ylim([0 1.5]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center_Relative.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(41);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin), Control');
xlim([0 60])
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_Cont_center.eps', '-depsc2');

figure(42);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, Control');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_Cont_center.eps', '-depsc2');

figure(43);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, Control');
xlim([0 60]);
ylim([0 1.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Cont_center.eps', '-depsc2');

figure(44);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, Relative Depth, Control');
xlim([0 5])
ylim([0 1.5]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Control_center_Relative.eps', '-depsc2');
%Control ends%
%}

%% Add labels (center slice) (Population mode) 

%GdNP+RT starts%
figure(11);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin, pop), GdNP+RT');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_GdNP+RT_center_pop.eps', '-depsc2');

figure(12);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, pop, GdNP+RT');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP+RT_center_pop.eps', '-depsc2');

figure(13);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, pop, GdNP+RT');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');

figure(14);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, pop, Relative Depth, GdNP+RT');
xlim([0 5])
ylim([0 0.35]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_center_pop_Relative.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(21);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin, pop), RT');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_RT_center_pop.eps', '-depsc2');

figure(22);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, pop, RT');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_RT_center_pop.eps', '-depsc2');

figure(23);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, pop, RT');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center_pop.eps', '-depsc2');

figure(24);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, pop, Relative Depth, RT');
xlim([0 5])
ylim([0 0.35]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center_pop_Relative.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(31);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin, pop), GdNP');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_GdNP_center_pop.eps', '-depsc2');

figure(32);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, pop, GdNP');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP_center_pop.eps', '-depsc2');

figure(33);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, pop, GdNP');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center_pop.eps', '-depsc2');

figure(34);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, pop, Relative Depth, GdNP');
xlim([0 5])
ylim([0 0.35]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center_pop_Relative.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(41);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans} (Large bin, pop), Control');
xlim([0 60])
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_LargeBin_Cont_center_pop.eps', '-depsc2');

figure(42);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Voxel count, Control, pop');
xlim([0 60])
ylim([0 300]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_Cont_center_pop.eps', '-depsc2');

figure(43);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('Vol-normalized K^{trans}, Control, pop');
xlim([0 60]);
ylim([0 0.35]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Cont_center_pop.eps', '-depsc2');

figure(44);
set(gca, 'LineWidth', 2, 'FontSize', 16);
title('Vol-normalized K^{trans}, pop, Relative Depth, Control');
xlim([0 5])
ylim([0 0.35]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Control_center_pop_Relative.eps', '-depsc2');
%Control ends%
