clear;
clc;
%% Prepare figure
figure(1); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(2); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(3); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(4); clf('reset');set(gcf,'Position',[440   378   560   420]);
figure(5); clf('reset');set(gcf,'Position',[440   378   560   420]);

%% 
%Convert voxel to volume
VoxelVolume = 0.15626 * 0.15625 * 1; % [mm]
basefolder = '/Users/Kwon/Documents/MATLAB/DCE_MRI/Analysis/Input';

AnimalID = [334 337 338 340 342 343 344 345 350 351 352 353 355];
AnimalID = AnimalID(:);
%1:334, 2:337, 3:338, 4:340, 5:342, 6:343, 7:344, 8:345, 9:350,
%10:351, 11:352, 12:353, 13:355
%GdNP+RT:1,2,4,8; RT:10,11,12; GdNP:3,6,7; Control:5,9,13

KtransFrac = zeros(13,100);
KtransNum = zeros(13,100);
TotalVoxelNum = zeros(13,1);
TotalKtrans = zeros(13,1);
AveKtrans = zeros(13,1);
AveKtransFracGdNPRT = zeros(13,100);
AveKtransFracRT = zeros(13,100);
AveKtransFracGdNP = zeros(13,100);
AveKtransFracControl = zeros(13,100);
%% Main Program
for k = 1:size(AnimalID,1)
    %Import Tumor mask
    fname_mask = sprintf('%d - as a 90 frames MultiVolume by ImagePositionPatientInstanceNumber frame 89-label.nrrd', AnimalID(k));
    fname = fullfile(basefolder, fname_mask);
    if exist(fname_mask,'file')
        [Vmask,info_mask] = nrrdread2(fname);
        %imshow(Vmask(:,:,2),[]);
        
        TotalVoxelNum(k,1) = sum(Vmask(:)); %Calculate total voxel number
        
        %Import Ktrans map
        fktrans = sprintf('%d_Output Ktrans image_pop.nrrd', AnimalID(k));
        fname_ktans = fullfile(basefolder, fktrans);
        if exist(fname_ktans,'file')
            [Vktrans,info_ktrans] = nrrdread2(fname_ktans);
            
            Bin = 5;
            for z = 1:5
                for y = 1:192
                    for x = 1:192
                        for Ktrans = 1:1:5*Bin
                            if((Ktrans-1)/Bin<Vktrans(x,y,z) && Vktrans(x,y,z)<=Ktrans/Bin)
                                KtransNum(k,Ktrans) = KtransNum(k,Ktrans) + 1;
                            end
                        end
                    end
                end
            end
            
            TotalKtrans(k,1) = sum(Vktrans(:));
            AveKtrans(k,1) = TotalKtrans(k,1)/TotalVoxelNum(k,1);
            
            %Convert number into fraction
            for Ktrans = 1:1:5*Bin
                KtransFrac(k,Ktrans) = KtransNum(k,Ktrans)/TotalVoxelNum(k,1);
            end
        else
            fprintf('Animal_%d, Ktrans map not imported\n', AnimalID(k));
        end        
    else
        fprintf('Animal_%d, mask not imported\n', AnimalID(k));
    end
end

for Ktrans = 1:1:5*Bin
    AveKtransFracGdNPRT(1,Ktrans) = (KtransFrac(1,Ktrans)+KtransFrac(2,Ktrans)+KtransFrac(4,Ktrans)+KtransFrac(8,Ktrans))/4;
    AveKtransFracRT(1,Ktrans) = (KtransFrac(10,Ktrans)+KtransFrac(11,Ktrans)+KtransFrac(12,Ktrans))/3;
    AveKtransFracGdNP(1,Ktrans) = (KtransFrac(3,Ktrans)+KtransFrac(6,Ktrans)+KtransFrac(7,Ktrans))/3;
    AveKtransFracControl(1,Ktrans) = (KtransFrac(5,Ktrans)+KtransFrac(9,Ktrans)+KtransFrac(13,Ktrans))/3;
end
%% Average Ktrans for entire volume
fprintf('Average Ktrans for GdNP+RT = %d\n', AveKtrans(1,1)+AveKtrans(2,1)+AveKtrans(4,1)+AveKtrans(8,1));
fprintf('Average Ktrans for RT = %d\n', AveKtrans(10,1)+AveKtrans(11,1)+AveKtrans(12,1));
fprintf('Average Ktrans for GdNP = %d\n', AveKtrans(3,1)+AveKtrans(6,1)+AveKtrans(7,1));
fprintf('Average Ktrans for Control = %d\n', AveKtrans(5,1)+AveKtrans(9,1)+AveKtrans(13,1));


%% Draw Figure
for k = 1:size(AnimalID,1)
    if (k == 1 || k == 2 || k == 4 || k == 8) %GdNP+RT
        if k == 1
            figure(1); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 2
            figure(1); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 4
            figure(1); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 8
            figure(1); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-^','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 10 || k == 11 || k == 12) %RT
        if k == 10
            figure(2); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 11
            figure(2); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 12
            figure(2); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 3 || k == 6 || k == 7) %GdNP
        if k == 3
            figure(3); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 6
            figure(3); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 7
            figure(3); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end

    elseif (k == 5 || k == 9 || k == 13) %Control
        if k == 5
            figure(4); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-o','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 9
            figure(4); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-+','LineWidth',2,'MarkerSize',10);
            hold off;
        elseif k == 13
            figure(4); hold on;
            plot(1/Bin:1/Bin:5,KtransFrac(k,1:5*Bin),'-*','LineWidth',2,'MarkerSize',10);
            hold off;
        end
    else
        fprintf('Animal_%d, mask not imported\n', AnimalID(k));
    end
end
figure(5); hold on;
plot(1/Bin:1/Bin:5,AveKtransFracControl(1,1:5*Bin),'-o','LineWidth',2,'MarkerSize',10); hold on;
plot(1/Bin:1/Bin:5,AveKtransFracGdNP(1,1:5*Bin),'-+','LineWidth',2,'MarkerSize',10); hold on;
plot(1/Bin:1/Bin:5,AveKtransFracRT(1,1:5*Bin),'-*','LineWidth',2,'MarkerSize',10); hold on;
plot(1/Bin:1/Bin:5,AveKtransFracGdNPRT(1,1:5*Bin),'-^','LineWidth',2,'MarkerSize',10); hold off;
%{
%% Add labels (AverageUnderAIFMask)
%GdNP+RT starts%
figure(1);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, GdNP+RT');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_GdNP+RT.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(2);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, RT');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_RT.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(3);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, GdNP');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_GdNP.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(4);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, Control');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_Control.eps', '-depsc2');
%Control ends%

%Average starts%
figure(5);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, Average');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Control','GdNP','RT','GdNP+RT','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_Ave.eps', '-depsc2');
%Control ends%
%}

%% Add labels (Population)
%GdNP+RT starts%
figure(1);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, pop, GdNP+RT');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal334','Animal337','Animal340','Animal345','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_GdNP+RT_pop.eps', '-depsc2');
%GdNP+RT ends%

%RT starts%
figure(2);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, pop, RT');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal351','Animal352','Animal353','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_RT_pop.eps', '-depsc2');
%RT ends%

%GdNP starts%
figure(3);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, pop, GdNP');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal338','Animal343','Animal344','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_GdNP_pop.eps', '-depsc2');
%GdNP ends%

%Control starts%
figure(4);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, pop, Control');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Animal342','Animal350','Animal355','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_Control_pop.eps', '-depsc2');
%Control ends%

%Average starts%
figure(5);
set(gca, 'LineWidth', 2, 'FontSize', 20);
title('K^{trans} fraction, pop, Average');
xlim([0 5]);
ylim([0 1]);
xlabel('K^{trans} (min^{-1})');
ylabel('Fraction');
legend('Control','GdNP','RT','GdNP+RT','Location','northeast');
legend boxoff;
% Save as eps
set(gcf, 'PaperPositionMode', 'Auto');
print(gcf, 'FracKtrans_Ave_pop.eps', '-depsc2');
%Control ends%