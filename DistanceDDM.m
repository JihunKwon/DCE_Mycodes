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
figure(51); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(52); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(53); clf('reset');set(gcf,'Position',[440   378   560   420]);

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
TotCount = zeros(13,1,DepthBinLimit);
Total = zeros(13,1,DepthBinLimit);
NetCount = zeros(13,1,DepthBinLimit);
Net = zeros(13,1,DepthBinLimit);
Ave = zeros(13,5,DepthBinLimit);

NetCountBin = zeros(13,1,20);
NetBin = zeros(13,1,20);
AveBin = zeros(13,1,20);
NetCountBinMax = zeros(13,20);
NetBinMax = zeros(13,20);
AveBinMax = zeros(13,20);

NetCountBinSiju = zeros(13,1,20);
NetBinSiju = zeros(13,1,20);
AveBinSiju = zeros(13,1,20);
NetCountBinMaxSiju = zeros(13,20);
NetBinMaxSiju = zeros(13,20);
AveBinMaxSiju = zeros(13,20);

TargetLog = zeros(13,1);
KtransFrac = zeros(13,DepthBinLimit);
AveKtrans = zeros(13,1);

NetCountBinRelative = zeros(13,1,20);
NetBinRelative = zeros(13,1,20);
AveBinRelative = zeros(13,1,20);
NetCountBinMaxRelative = zeros(13,20);
NetBinMaxRelative = zeros(13,20);
AveBinMaxRelative = zeros(13,20);

MaxDepth5 = zeros(13,6); %Relative depth from edge. (AnimalID, Slice, MaxDepth)
RelativeCount5 = zeros(13,10);
RelativeNet5 = zeros(13,10);
RelativeAve5 = zeros(13,10);

MaxDepth10 = zeros(13,6); %Relative depth from edge. (AnimalID, Slice, MaxDepth)
RelativeCount10 = zeros(13,10);
RelativeNet10 = zeros(13,10);
RelativeAve10 = zeros(13,10);

AveAveBinMaxSiju = zeros(10,4);
AveAveBinMax = zeros(10,4);
AveRelativeAve = zeros(10,4);

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
            fktrans = sprintf('%d_Output Ktrans image.nrrd', AnimalID(k));
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
                %For 20% bin 
                for i=DepthBinLimit:-1:1
                    depth = FindMaxDepth(Net,k,i); %Find Maximum i. i=depth
                    if depth~=0
                        MaxDepth5(k,1) = depth;
                        break;
                    end
                end
                
                for j = 2:6
                    MaxDepth5(k,j) = MaxDepth5(k,1)/100*(j-1)*20; %Convert to 20%,40%,...,100%
                end
                MaxDepth5(k,1)=0; %To use this variable in the next for loop
                
                for i = 1:1:DepthBinLimit
                    if i <= MaxDepth5(k,2)
                        RelativeCount5(k,1) = RelativeCount5(k,1) + NetCount(k,1,i);
                        RelativeNet5(k,1) = RelativeNet5(k,1) + Net(k,1,i);
                    elseif i <= MaxDepth5(k,3)
                        RelativeCount5(k,2) = RelativeCount5(k,2) + NetCount(k,1,i);
                        RelativeNet5(k,2) = RelativeNet5(k,2) + Net(k,1,i);
                    elseif i <= MaxDepth5(k,4)
                        RelativeCount5(k,3) = RelativeCount5(k,3) + NetCount(k,1,i);
                        RelativeNet5(k,3) = RelativeNet5(k,3) + Net(k,1,i);
                    elseif i <= MaxDepth5(k,5)
                        RelativeCount5(k,4) = RelativeCount5(k,4) + NetCount(k,1,i);
                        RelativeNet5(k,4) = RelativeNet5(k,4) + Net(k,1,i);
                    elseif i <= MaxDepth5(k,6)
                        RelativeCount5(k,5) = RelativeCount5(k,5) + NetCount(k,1,i);
                        RelativeNet5(k,5) = RelativeNet5(k,5) + Net(k,1,i);
                    end
                end
                
                %For 10% bin
                for i=DepthBinLimit:-1:1
                    depth = FindMaxDepth(Net,k,i); %Find Maximum i. i=depth
                    if depth~=0
                        MaxDepth10(k,1) = depth;
                        break;
                    end
                end
                for j = 2:11
                    MaxDepth10(k,j) = MaxDepth10(k,1)/100*(j-1)*10; %Convert to 10%,20%,...,100%
                end
                MaxDepth10(k,1)=0; %To use this variable in the next for loop
                
                for i = 1:1:DepthBinLimit
                    if i <= MaxDepth10(k,2)
                        RelativeCount10(k,1) = RelativeCount10(k,1) + NetCount(k,1,i);
                        RelativeNet10(k,1) = RelativeNet10(k,1) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,3)
                        RelativeCount10(k,2) = RelativeCount10(k,2) + NetCount(k,1,i);
                        RelativeNet10(k,2) = RelativeNet10(k,2) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,4)
                        RelativeCount10(k,3) = RelativeCount10(k,3) + NetCount(k,1,i);
                        RelativeNet10(k,3) = RelativeNet10(k,3) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,5)
                        RelativeCount10(k,4) = RelativeCount10(k,4) + NetCount(k,1,i);
                        RelativeNet10(k,4) = RelativeNet10(k,4) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,6)
                        RelativeCount10(k,5) = RelativeCount10(k,5) + NetCount(k,1,i);
                        RelativeNet10(k,5) = RelativeNet10(k,5) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,7)
                        RelativeCount10(k,6) = RelativeCount10(k,6) + NetCount(k,1,i);
                        RelativeNet10(k,6) = RelativeNet10(k,6) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,8)
                        RelativeCount10(k,7) = RelativeCount10(k,7) + NetCount(k,1,i);
                        RelativeNet10(k,7) = RelativeNet10(k,7) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,9)
                        RelativeCount10(k,8) = RelativeCount10(k,8) + NetCount(k,1,i);
                        RelativeNet10(k,8) = RelativeNet10(k,8) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,10)
                        RelativeCount10(k,9) = RelativeCount10(k,9) + NetCount(k,1,i);
                        RelativeNet10(k,9) = RelativeNet10(k,9) + Net(k,1,i);
                    elseif i <= MaxDepth10(k,11)
                        RelativeCount10(k,10) = RelativeCount10(k,10) + NetCount(k,1,i);
                        RelativeNet10(k,10) = RelativeNet10(k,10) + Net(k,1,i);
                    end
                end
                
                for i = 1:1:DepthBinLimit
                    for j = 1:5
                        RelativeAve5(k,j) = RelativeNet5(k,j)/RelativeCount5(k,j);
                    end
                    for j = 1:10
                        RelativeAve10(k,j) = RelativeNet10(k,j)/RelativeCount10(k,j);
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

%1:NP+RT, 2:RT, 3:Gd, 4:Control
AveAveBinMaxSiju(1,1:RoundMaxSiju) = (AveBinMaxSiju(1,1:RoundMaxSiju)+AveBinMaxSiju(2,1:RoundMaxSiju)+AveBinMaxSiju(3,1:RoundMaxSiju)+AveBinMaxSiju(8,1:RoundMaxSiju))/4;
AveAveBinMaxSiju(2,1:RoundMaxSiju) = (AveBinMaxSiju(10,1:RoundMaxSiju)+AveBinMaxSiju(11,1:RoundMaxSiju)+AveBinMaxSiju(12,1:RoundMaxSiju))/3;
AveAveBinMaxSiju(3,1:RoundMaxSiju) = (AveBinMaxSiju(3,1:RoundMaxSiju)+AveBinMaxSiju(6,1:RoundMaxSiju)+AveBinMaxSiju(7,1:RoundMaxSiju))/3;
AveAveBinMaxSiju(4,1:RoundMaxSiju) = (AveBinMaxSiju(5,1:RoundMaxSiju)+AveBinMaxSiju(9,1:RoundMaxSiju)+AveBinMaxSiju(13,1:RoundMaxSiju))/3;

AveAveBinMax(1,1:RoundMax) = (AveBinMax(1,1:RoundMax)+AveBinMax(2,1:RoundMax)+AveBinMax(3,1:RoundMax)+AveBinMax(8,1:RoundMax))/4;
AveAveBinMax(2,1:RoundMax) = (AveBinMax(10,1:RoundMax)+AveBinMax(11,1:RoundMax)+AveBinMax(12,1:RoundMax))/3;
AveAveBinMax(3,1:RoundMax) = (AveBinMax(3,1:RoundMax)+AveBinMax(6,1:RoundMax)+AveBinMax(7,1:RoundMax))/3;
AveAveBinMax(4,1:RoundMax) = (AveBinMax(5,1:RoundMax)+AveBinMax(9,1:RoundMax)+AveBinMaxSiju(13,1:RoundMax))/3;

AveRelativeAve5(1,1:5) = (RelativeAve5(1,1:5)+RelativeAve5(2,1:5)+RelativeAve5(3,1:5)+RelativeAve5(8,1:5))/4;
AveRelativeAve5(2,1:5) = (RelativeAve5(10,1:5)+RelativeAve5(11,1:5)+RelativeAve5(12,1:5))/3;
AveRelativeAve5(3,1:5) = (RelativeAve5(3,1:5)+RelativeAve5(6,1:5)+RelativeAve5(7,1:5))/3;
AveRelativeAve5(4,1:5) = (RelativeAve5(5,1:5)+RelativeAve5(9,1:5)+RelativeAve5(13,1:5))/3;

AveRelativeAve10(1,1:10) = (RelativeAve10(1,1:10)+RelativeAve10(2,1:10)+RelativeAve10(3,1:10)+RelativeAve10(8,1:10))/4;
AveRelativeAve10(2,1:10) = (RelativeAve10(10,1:10)+RelativeAve10(11,1:10)+RelativeAve10(12,1:10))/3;
AveRelativeAve10(3,1:10) = (RelativeAve10(3,1:10)+RelativeAve10(6,1:10)+RelativeAve10(7,1:10))/3;
AveRelativeAve10(4,1:10) = (RelativeAve10(5,1:10)+RelativeAve10(9,1:10)+RelativeAve10(13,1:10))/3;

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
            plot(1:10,RelativeAve10(k,1:10),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 2
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 4
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 8
            figure(11); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-^','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(1,1:RoundMaxSiju),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(12); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(13); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-^','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:10,AveRelativeAve10(1,1:10),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(14); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-^','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:5,AveRelativeAve5(1,1:5),'--kd','LineWidth',2,'MarkerSize',10);
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
            plot(1:10,RelativeAve10(k,1:10),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 11
            figure(21); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(22); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(23); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 12
            figure(21); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(2,1:RoundMaxSiju),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(22); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(23); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:10,AveRelativeAve10(2,1:10),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(24); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:5,AveRelativeAve5(2,1:5),'--kd','LineWidth',2,'MarkerSize',10);
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
            plot(1:10,RelativeAve10(k,1:10),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 6
            figure(31); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(32); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(33); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 7
            figure(31); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(3,1:RoundMaxSiju),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(32); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(33); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:10,AveRelativeAve10(3,1:10),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(34); hold on;
            plot(1:5,RelativeAve5(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:5,AveRelativeAve5(3,1:5),'--kd','LineWidth',2,'MarkerSize',10);
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
            plot(1:10,RelativeAve10(k,1:10),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 9
            figure(41); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(42); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(43); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 13
            figure(41); hold on;
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveBinMaxSiju(k,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(4,1:RoundMaxSiju),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(42); hold on;
            plot(ModulateBin:ModulateBin:MaxBin,NetCountBinMax(k,1:RoundMax),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(43); hold on;
            plot(1:10,RelativeAve10(k,1:10),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:10,AveRelativeAve10(4,1:10),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
            figure(44); hold on;
            plot(1:1:5,RelativeAve5(k,1:5),'-*','LineWidth',2,'MarkerSize',10);
            hold on; %Average
            plot(1:5,AveRelativeAve5(4,1:5),'--kd','LineWidth',2,'MarkerSize',10);
            hold off;
        end
    else
        fprintf('Animal_%d, mask not imported\n', AnimalID(k));
    end
end

%plot Averages
figure(51); hold on;
plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(4,1:RoundMaxSiju),'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(3,1:RoundMaxSiju),'-+','LineWidth',2,'MarkerSize',10); hold on;
plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(2,1:RoundMaxSiju),'-*','LineWidth',2,'MarkerSize',10); hold on;
plot(ModulateBinSiju:ModulateBinSiju:MaxBin,AveAveBinMaxSiju(1,1:RoundMaxSiju),'-^','LineWidth',2,'MarkerSize',10); hold on;
hold off;
figure(52); hold on;
plot(1:10,AveRelativeAve10(4,1:10),'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(1:10,AveRelativeAve10(3,1:10),'-+','LineWidth',2,'MarkerSize',10); hold on;
plot(1:10,AveRelativeAve10(2,1:10),'-*','LineWidth',2,'MarkerSize',10); hold on;
plot(1:10,AveRelativeAve10(1,1:10),'-^','LineWidth',2,'MarkerSize',10); hold on;
hold off;
figure(53); hold on;
plot(1:5,AveRelativeAve5(4,1:5),'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(1:5,AveRelativeAve5(3,1:5),'-+','LineWidth',2,'MarkerSize',10); hold on;
plot(1:5,AveRelativeAve5(2,1:5),'-*','LineWidth',2,'MarkerSize',10); hold on;
plot(1:5,AveRelativeAve5(1,1:5),'-^','LineWidth',2,'MarkerSize',10); hold on;
hold off;
%% Add labels (center slice) (Population mode) 

%GdNP+RT starts%
figure(11);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP+RT');
xlim([0 60])
ylim([0 0.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_LargeBin_center.eps', '-depsc2');

figure(12);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Voxel count, GdNP+RT');
xlim([0 60])
ylim([0 500]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP+RT_center.eps', '-depsc2');

figure(13);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP+RT');
xlim([1 10]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_center_Relative_10%.eps', '-depsc2');

figure(14);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP+RT');
xlim([1 5]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP+RT_center_Relative_20%.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(21);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, RT');
xlim([0 60])
ylim([0 0.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_LargeBin_center.eps', '-depsc2');

figure(22);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Voxel count, RT');
xlim([0 60])
ylim([0 400]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_RT_center.eps', '-depsc2');

figure(23);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, RT');
xlim([1 10]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center_Relative_10%.eps', '-depsc2');

figure(24);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, RT');
xlim([1 5]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_RT_center_Relative_20%.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(31);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP');
xlim([0 60])
ylim([0 0.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_LargeBin_center.eps', '-depsc2');

figure(32);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Voxel count, GdNP');
xlim([0 60])
ylim([0 400]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_GdNP_center.eps', '-depsc2');

figure(33);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP');
xlim([1 10]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center_Relative_10%.eps', '-depsc2');

figure(34);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, GdNP');
xlim([1 5]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_GdNP_center_Relative_20%.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(41);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Control');
xlim([0 60])
ylim([0 0.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Cont_LargeBin_center.eps', '-depsc2');

figure(42);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Voxel count, Control');
xlim([0 60])
ylim([0 400]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Voxel count');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'Count_Cont_center.eps', '-depsc2');

figure(43);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Control');
xlim([1 10]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Control_center_Relative_10%.eps', '-depsc2');

figure(44);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Control');
xlim([1 5]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Control_center_Relative_20%.eps', '-depsc2');
%Control ends%

%Average starts%
figure(51);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Average');
xlim([0 60])
ylim([0 0.5]);
xlabel('Depth from the tumor surface (mm)');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('GdNP+RT','RT-only','GdNP-only','Control','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Ave_LargeBin_center.eps', '-depsc2');

figure(52);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Average');
xlim([0 10]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Control','GdNP-only','RT-only','GdNP+RT','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Ave_center_Relative_10%.eps', '-depsc2');
%Average ends%

figure(53);
set(gca, 'LineWidth', 2, 'FontSize', 18);
title('Vol-normalized K^{trans}, Average');
xlim([0 5]);
ylim([0 0.4]);
xlabel('Relative distance from tumor edge');
ylabel('Vol-normalized K^{trans} (min^{-1})');
legend('Control','GdNP-only','RT-only','GdNP+RT','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'KtransAve_Ave_center_Relative_20%.eps', '-depsc2');
%Average ends%
