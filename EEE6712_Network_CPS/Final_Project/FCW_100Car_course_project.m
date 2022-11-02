% Copyright 2017 ; Yaser P. Fallah, UCF, Dept. of Electrical and Computer Engineering

function [] = FCW_100Car_course_project_v7()
clear all; close all; clc;
timestep = 0.1; %sec
G = 9.8;
DEBUG = 0;
DO_EDCFS = 0;
filenumbers = 8296:1:9122;
global CONST_SPEED_MODEL;  % 1 means constat speed model, 0 means const. acceleration model.
CONST_SPEED_MODEL = 0;

%--------read from files, and process for alert
%PAPERS:   look at "Vehicle lateral and longitudinal velocity estimation based on Unscented Kalman Filter"

% fnameprefix = '..\\100CarData\\100CarAllData_v1_3\\HundredCar_Public_%d.txt';
 fnameprefix = '/home/nitish/Documents/FCWProject/100CarData/100CarAllData_v1_3/HundredCar_Public_%d.txt';
% fnameprefix = '/home/nitish/Desktop/UCFCourses/Spring2019/EEE6712/FCWProject/FCWProject/100CarData/100CarAllData_v1_3/HundredCar_Public_%d.txt';

PER_values = zeros(1,10);
PER_counter = 0;
%meanAcc = 0;
for PER =0:0.1:0.9
    PER_counter = PER_counter+1;
    Rate = 10;
    PER_values(PER_counter)=PER;
    scenario = 0;
    AnalysisResult_PB_camp = zeros(1,7);
    AnalysisResult_PB_knipling = zeros(1,7);
    AnalysisResult_PB_NHTSA_early = zeros(1,7);
    AnalysisResult_PB_NHTSA_inter = zeros(1,7);
    AnalysisResult_PB_NHTSA_immi = zeros(1,7);
    
    for filenum = 8300:1:8400;%filenumbers
    FileNumber = filenum
    scenario = scenario + 1;
    fname = sprintf(fnameprefix,filenum);
    [fid,msg] = fopen(fname);
    C = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %*[^\n]', 'delimiter', ',');
    fclose(fid); tripid = C{1}; sync = C{2}; triptime = C{3}; gaspedpos=C{4}; speedvc=C{5};
    speedgpshrz=C{6}; yaw=C{7} ; gpsheading = C{8}; latAcc = C{9}; longAcc = C{10};
    %if(DEBUG) figure();plot(speedvc);end;
    
    %initialize LV and FV
    N=size(triptime, 1);
    x_lv=zeros(1,N);v_lv=zeros(1,N);accel= zeros(1,N);a_lv_calc = zeros(1,N);
    x_fv=zeros(1,N);v_fv=zeros(1,N);a_fv = zeros(1,N);
    Rw=zeros(1,N); %warning range
    GrndTrth_Separation_Dist = zeros(1,N);
    RwCamp = zeros(1,N);
    RwCamph = zeros(1,N);
    %RwKnipling = zeros(1,N);
    
    v_lv = 0.44704*(speedvc'); %change to m/s
    if (v_lv(1) < 0) v_lv(1) = 0; end;
    
    % Cleaning acceleration data, removing bias and smoothing acceleration
    acclong = G*(longAcc);% - mean(longAcc(1:69)) );
    acclat = G*(latAcc);%  - mean(latAcc(1:25)) );
    heading = 0.1*(yaw - mean(yaw(1:25)));
    accel = acclong.*cos(heading)-acclat.*sin(heading);
    nn = min (N,150);
    acc_bias = (   sum(accel(1:nn)) *timestep - (v_lv(nn)-v_lv(1)) )/(nn*timestep);
    a_lv_unbiased = accel - acc_bias;
    a_lv = filterAcc(a_lv_unbiased, 5); % odd number averaging filter
    %      figure(45); hold on; plot(a_lv);plot(a_lv_unbiased,'r');
    
    v_lv_calc(1) = v_lv(1);
    %v_lv_calc
    tripsamples = size (triptime,1);
    v_fv(1) = v_lv(1) ;%m/s
    x_fv(1) = 0; x_lv(1)= x_fv(1) + 40;
    cfmode = zeros(1,N);
    
    % Car following model from MITSIM is implemented here. It calculates
    % following vehicle (host vehicle) position and speed, using
    % actual information from lead vehicle (remote vehicle)
    for t=2:tripsamples
        if v_lv(t) < 0   % for the values of V that are missing, -1 is entered in VTTI data, use the previous speed instead
            v_lv(t) = v_lv(t-1);
        end
        % should we use constant acceleraton or constant speed model for
        a_lv_calc(t) = (v_lv(t) - v_lv(t-1))/timestep ;
        v_lv_calc(t) = v_lv_calc(t-1)+ a_lv(t-1)*timestep + 0.02*(0.5-rand(1));
        v_lv(t) = v_lv_calc(t) ;
        x_lv(t) = x_lv(t-1) + v_lv(t-1)*timestep + 0.5*a_lv(t-1)*timestep^2 + 0.02*(0.5-rand(1)); %assume constant acceleration from previous sample + 0.5*a_lv(t-1)*timestep^2; %assume constant acceleration from previous sample
        %v_lv(t) = v_lv(t-1) + a_lv(t-1)*timestep;
        
        x_fv(t) = x_fv(t-1) + v_fv(t-1)*timestep + 0.5*a_fv(t-1)*timestep^2; %assume constant acceleration from previous sample
        v_fv(t) = v_fv(t-1) + a_fv(t-1)*timestep;
        
        if ( v_fv(t) <= 0.44704*0)
            v_fv(t) = 0.44704*0.0;
            %         'FV stopped at '
            %         t
        end
        
        %car-following model
        VehLen = 4.5; %meters
        HlowerT = 0.5; HupperT = 1.36; Vdesired = 0.44704*65;
        Hupper = HupperT * ( v_fv(t) ) + 20 ;
        Hlower = HlowerT * ( v_fv(t) ) + 5;
        if( x_lv(t) - x_fv(t) -VehLen > Hupper )  %free driving
            cfmode(t) = 1;
            if (v_fv(t) < Vdesired)
                a_fv(t) = aplusf(v_fv(t));
            elseif v_fv(t) == Vdesired
                a_fv(t) = 0;
            else
                a_fv(t) = aminusf(v_fv(t));
            end
            
        elseif( x_lv(t) - x_fv(t) -VehLen < Hlower )  %emergency
            decel = 1*G; %.5G
            cfmode(t) = -1;
            if( x_lv(t) - x_fv(t) -VehLen <= 0 )
                if (DEBUG)Hlower,Hupper, end;
                separation_distance = x_lv(t) - x_fv(t)  -VehLen;
                x_fv(t) = x_lv(t) - VehLen -1;
            else
                if( v_fv(t) > v_lv(t) )
                    decel = min(aminusf(v_fv(t)), (a_lv(t)- 0.5*((v_fv(t) - v_lv(t))^2) / (x_lv(t) - VehLen - x_fv(t)) ) );
                else
                    decel = min(aminusf(v_fv(t)), a_lv(t) + 0.25*aminusf(v_fv(t)));
                end
            end
            a_fv(t) = decel;
            
        else %car following if ( x_lv(t) - x_fv(t) > 100 )
            cfmode(t) = 0;
            if ( v_fv(t) > v_lv(t) ) %deceleration
                alpha = 1.55;beta = 1.08; gamma = 1.65;
            else %acceleration
                alpha = 2.15;beta = -1.67; gamma = -0.89;
            end
            if( v_fv(t) > 0.44*2 )
                %Yang's formula
                a_fv(t) = alpha*(v_fv(t)^beta)*( (v_lv(t) - v_fv(t)) / ((x_lv(t) - VehLen- x_fv(t))^gamma) );
                if(a_fv(t)>0.5*G )
                    a_fv(t) = 0.5*G;
                elseif (a_fv(t) < -0.8*G)
                    a_fv(t) = -0.8*G;
                end
            else
                a_fv(t) = a_lv(t);
            end
            
        end
        GrndTrth_Separation_Dist(t) = x_lv(t) - x_fv(t) - VehLen;
        
        
        
    end%end of trip,
    
         %figure(); hold on; plot(a_lv, 'b');plot(a_lv_calc,'r');plot(acclong,'k'); plot(acclat,'g'); %plot(a_lv,'b');
         %figure();plot(v_lv); hold on; plot(v_lv_calc,'r');

    triptime = 0:timestep:(tripsamples*timestep)-timestep;
    
    %Warning Algorithm, establish ground truth for ideal communication case
    RwCamp = CAMPLinearWarningAlg(v_lv,v_fv,a_lv,a_fv); %returns warning range
    GroundTruthAlertRange = RwCamp;
    %plot(triptime,Range,'b');
    %figure(1);plot(triptime,RwCamp,'r');xlabel('triptime');ylabel('Rw');hold on;

    
    % sample test with a particular Rate and PER value. This function
    % returns the received information of LV (remote vehicle) over a
    % network (using periodic sampling pf LV position). PBRate will be
    % almost Rate*(1-PER), 
    %Rate = 10; PER = 0.90; %at 82%, camp is generating more warnings than knipling
    %PER = 0.30;
    [v_lv_rcvd, a_lv_rcvd, x_lv_rcvd, PBRate] = findReceivedLVInfoOverPBNetwork(v_lv,a_lv,x_lv,Rate,PER);
    
    %position tracking error - error in position of LV as calculated in the
    %FV (host vehicle)
    x_lv_est =  x_lv_rcvd;
    PTE(scenario)=prctile (abs (x_lv - x_lv_est),95);
    est_Separation_Dist = x_lv_est - x_fv - VehLen;
    
    WarningRangeCalculatedUsingReceivedInfo = CAMPLinearWarningAlg(v_lv_rcvd,v_fv,a_lv_rcvd, a_fv);
    [RwKnipling,lvstate] = Knipling(v_lv_rcvd,v_fv,a_lv_rcvd, a_fv);
    %maxrv
    %maxVfv
    %maxAfv
    %maxVlv
    %maxAlv
    %max_val = max(RwKnipling)

    
    %[Dmiss_early_out,Dmiss_inter_out,Dmiss_immi_out] = NHTSA(x_lv_est,x_fv,VehLen,v_lv_rcvd,v_fv,a_lv_rcvd, a_fv);
    [Dmiss_early_out,Dmiss_inter_out,Dmiss_immi_out] = NHTSA_2(x_lv_est,x_fv,VehLen,v_lv_rcvd,v_fv,a_lv_rcvd, a_fv);
%     early_warning_counter
%     inter_warning_counter
%     immi_warning_counter
%     no_warning_counter
    
    for i=1:size(v_fv,2)
        NHTSA_thresh(i)= 2+v_fv(i)*0.1;
    end

%      warning_range_plot = figure();
%      xlabel('triptime[s]');ylabel('Distance[m]');hold on; plot(triptime,RwCamp,'r'); plot(triptime,RwKnipling,'g'); plot(triptime, est_Separation_Dist, 'k');legend('CAMP linear','Knipling', 'distance b/w vehicles');
%      Ttitle='Warning range vs triptime';
%      path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
%      saveas(warning_range_plot,path);
% 
%      Dmiss_plot = figure();
%      xlabel('triptime[s]');ylabel('Distance[m]');hold on;plot(triptime,Dmiss_early_out,'m');plot(triptime,Dmiss_inter_out,'b');plot(triptime,Dmiss_immi_out,'g');plot(triptime, est_Separation_Dist, 'k');legend('NHTSA early','NHTSA inter','NHTSA immi','distance b/w vehicles');
%      Ttitle='NHTSA Distance-to-Miss vs triptime';
%      path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
%      saveas(Dmiss_plot,path);

    num_warnings_camp = 0;
    num_warnings_knipling = 0;
    num_warnings_NHTSA_early=0;
    num_warnings_NHTSA_inter=0;
    num_warnings_NHTSA_immi=0;

    threatvec_early = zeros(1,size(triptime,2));
    threatvec_inter = zeros(1,size(triptime,2));
    threatvec_immi = zeros(1,size(triptime,2));

    safevec_early = zeros(1,size(triptime,2));
    safevec_inter = zeros(1,size(triptime,2));
    safevec_immi = zeros(1,size(triptime,2));

    for tt= 1:size(triptime, 2)
        if(est_Separation_Dist(tt) <= WarningRangeCalculatedUsingReceivedInfo(tt))
            num_warnings_camp = num_warnings_camp+1;
        end

        if(lvstate(tt)==0)
            if((est_Separation_Dist(tt)<=RwKnipling(tt)) && (v_lv_rcvd(tt)<(1609.34/3600)) && (v_fv(tt)>v_lv_rcvd(tt)))
                num_warnings_knipling = num_warnings_knipling+1;
%             else
%                 RwKnipling(tt)=0;
            end
        elseif(lvstate(tt)==1)
            if((est_Separation_Dist(tt)<=RwKnipling(tt)) && (v_lv_rcvd(tt)<v_fv(tt)) && (a_lv_rcvd(tt)<0))
                num_warnings_knipling = num_warnings_knipling+1;
%             else
%                 RwKnipling(tt)=0;
            end
        end 

        if(Dmiss_early_out(tt) <= 2+v_fv(tt)*0.1)
            num_warnings_NHTSA_early = num_warnings_NHTSA_early+1;
            threatvec_early(tt)=1;
            safevec_early(tt)=0;
        else
            threatvec_early(tt)=0;
            safevec_early(tt)=1;
        end
        if(Dmiss_inter_out(tt) <= 2+v_fv(tt)*0.1)
            num_warnings_NHTSA_inter = num_warnings_NHTSA_inter+1;
            threatvec_inter(tt)=1;
            safevec_inter(tt)=0;
        else
            threatvec_inter(tt)=0;
            safevec_inter(tt)=1;
        end
        if(Dmiss_immi_out(tt) <= 2+v_fv(tt)*0.1)
            num_warnings_NHTSA_immi = num_warnings_NHTSA_immi+1;
            threatvec_immi(tt)=1;
            safevec_immi(tt)=0;
        else
            threatvec_immi(tt)=0;
            safevec_immi(tt)=1;
        end
    end
    
    AnalysisResult_warnings_camp(scenario) = num_warnings_camp;
    AnalysisResult_warnings_knipling(scenario) = num_warnings_knipling;
    AnalysisResult_warnings_NHTSA_early(scenario) = num_warnings_NHTSA_early;
    AnalysisResult_warnings_NHTSA_inter(scenario) = num_warnings_NHTSA_inter;
    AnalysisResult_warnings_NHTSA_immi(scenario) = num_warnings_NHTSA_immi;
%     num_warnings_camp
%     num_warnings_knipling
%     num_warnings_NHTSA_early
%     num_warnings_NHTSA_inter
%     num_warnings_NHTSA_immi
    
    [a, b, c, d, e, f, g] = analyzeaccuracy(WarningRangeCalculatedUsingReceivedInfo, GroundTruthAlertRange, est_Separation_Dist,GrndTrth_Separation_Dist);
    AnalysisResult_PB_camp(scenario,1:7) = [a b c d e f g];
    [a, b, c, d, e, f, g] = analyzeaccuracy(RwKnipling, GroundTruthAlertRange, est_Separation_Dist,GrndTrth_Separation_Dist);
    AnalysisResult_PB_knipling(scenario,1:7) = [a b c d e f g];
    [a, b, c, d, e, f, g] = analyzeaccuracy_NHTSA(threatvec_early, safevec_early,GroundTruthAlertRange, est_Separation_Dist,GrndTrth_Separation_Dist);
    AnalysisResult_PB_NHTSA_early(scenario,1:7) = [a b c d e f g];
    [a, b, c, d, e, f, g] = analyzeaccuracy_NHTSA(threatvec_inter, safevec_inter,GroundTruthAlertRange, est_Separation_Dist,GrndTrth_Separation_Dist);
    AnalysisResult_PB_NHTSA_inter(scenario,1:7) = [a b c d e f g];
    [a, b, c, d, e, f, g] = analyzeaccuracy_NHTSA(threatvec_immi, safevec_immi,GroundTruthAlertRange, est_Separation_Dist,GrndTrth_Separation_Dist);
    AnalysisResult_PB_NHTSA_immi(scenario,1:7) = [a b c d e f g];
    end

    %accuracy
    meanAcc_camp(PER_counter) = mean(AnalysisResult_PB_camp(:,1));
    meanAcc_knipling(PER_counter) = mean(AnalysisResult_PB_knipling(:,1));
    meanAcc_NHTSA_early(PER_counter)= mean(AnalysisResult_PB_NHTSA_early(:,1));
    meanAcc_NHTSA_inter(PER_counter)= mean(AnalysisResult_PB_NHTSA_inter(:,1));
    meanAcc_NHTSA_immi(PER_counter)= mean(AnalysisResult_PB_NHTSA_immi(:,1));

    %precision
    meanPrec_camp(PER_counter) = mean(AnalysisResult_PB_camp(:,2));
    meanPrec_knipling(PER_counter) = mean(AnalysisResult_PB_knipling(:,2));
    meanPrec_NHTSA_early(PER_counter)= mean(AnalysisResult_PB_NHTSA_early(:,2));
    meanPrec_NHTSA_inter(PER_counter)= mean(AnalysisResult_PB_NHTSA_inter(:,2));
    meanPrec_NHTSA_immi(PER_counter)= mean(AnalysisResult_PB_NHTSA_immi(:,2));

    %warnings
    meanNumWarnings_camp(PER_counter) = mean(AnalysisResult_warnings_camp);
    meanNumWarnings_knipling(PER_counter) = mean(AnalysisResult_warnings_knipling);
    meanNumWarnings_NHTSA_early(PER_counter) = mean(AnalysisResult_warnings_NHTSA_early);
    meanNumWarnings_NHTSA_inter(PER_counter) = mean(AnalysisResult_warnings_NHTSA_inter);
    meanNumWarnings_NHTSA_immi(PER_counter) = mean(AnalysisResult_warnings_NHTSA_immi);

%     figure();hold on;xlabel('rate of transmitted information (Hz)');ylabel('Accuracy');
%          meanGmeanForEachRateBL(1:size(Raterange'.*(1-PERrange)',1))=mean(AnalysisResult_BL(:,1:size(PERrange,2),1:size(Raterange,2),7));
%          semilogx( Raterange'.*(1-PERrange)' , meanGmeanForEachRateBL );
%     xlim([0.2 10]);
end

% h=figure('pos', [10 10 1000 1200]);
% fontSize = 12;
% linewidthvalue = 1.1;
% [ha, pos] = tight_subplot(3, 1, [0.07 0.07], [0.07 0.07], [0.09 0.05]);
% set(gca,'FontSize', fontSize);
%  axes(ha(1));
% plot(PER_values,meanAcc_camp,'r');legend('Camp_linear');
% hold on;
%   ax = gca;
%             ax.XTick=0:0.1:0.9;
%             ax.XTickLabel=0:0.1:0.9;
%             ax.YTick=0:0.1:0.9;
%             ax.YTickLabel=0:0.1:0.9;
% ylim([0 1]);
%                    xlim([0 1]);
% xlabel('PER');ylabel('Mean Accuracy');
%  axes(ha(2));
% plot(PER_values,meanAcc_knipling,'b');legend('Knipling');
%   ax = gca;
%            ax.XTick=0:0.1:0.9;
%             ax.XTickLabel=0:0.1:0.9;
%             ax.YTick=0:0.1:0.9;
%             ax.YTickLabel=0:0.1:0.9;
% ylim([0 1]);
%                    xlim([0 1]);
% xlabel('PER');ylabel('Mean Accuracy');
% axes(ha(3));
% plot(PER_values,meanAcc_NHTSA_early,'m');plot(PER_values,meanAcc_NHTSA_inter,'b');plot(PER_values,meanAcc_NHTSA_immi,'r');legend('NHSTA early','NHSTA inter','NHSTA immi');
%   ax = gca;
%              ax.XTick=0:0.1:0.9;
%             ax.XTickLabel=0:0.1:0.9;
%             ax.YTick=0:0.1:0.9;
%             ax.YTickLabel=0:0.1:0.9;
% ylim([0 1]);
%                    xlim([0 1]);
% xlabel('PER');ylabel('Mean Accuracy');
% set(gca,'FontSize', fontSize);

% Ttitle='test';
% path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
% saveas(h,path)      



% figure();xlabel('PER');ylabel('Mean Accuracy');hold on;plot(PER_values,meanAcc_camp,'r');legend('Camp_linear');     
% figure();xlabel('PER');ylabel('Mean Accuracy');hold on;plot(PER_values,meanAcc_knipling,'b');legend('Knipling');
% figure();xlabel('PER');ylabel('Mean Accuracy');hold on;plot(PER_values,meanAcc_NHTSA_early,'m');plot(PER_values,meanAcc_NHTSA_inter,'b');plot(PER_values,meanAcc_NHTSA_immi,'r');legend('NHSTA early','NHSTA inter','NHSTA immi');

   

accuracy_plot = figure();
xlabel('PER');ylabel('Mean Accuracy');hold on;plot(PER_values,meanAcc_camp,'r');legend('Camp_linear');     
hold on;plot(PER_values,meanAcc_knipling,'b');legend('Knipling');
hold on;plot(PER_values,meanAcc_NHTSA_early,'m');plot(PER_values,meanAcc_NHTSA_inter,'g');plot(PER_values,meanAcc_NHTSA_immi,'r');legend('Camp linear', 'Knipling','NHSTA early','NHSTA inter','NHSTA immi');                

Ttitle='Accuracy_vs_PER';
path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
saveas(accuracy_plot,path);


precision_plot = figure();
xlabel('PER');ylabel('Mean Precision');hold on;plot(PER_values,meanPrec_camp,'r');legend('Camp_linear');     
hold on;plot(PER_values,meanPrec_knipling,'b');legend('Knipling');
hold on;plot(PER_values,meanPrec_NHTSA_early,'m');plot(PER_values,meanPrec_NHTSA_inter,'g');plot(PER_values,meanPrec_NHTSA_immi,'r');legend('Camp linear', 'Knipling','NHSTA early','NHSTA inter','NHSTA immi');                

Ttitle='Precision_vs_PER';
path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
saveas(precision_plot,path);


warnings_plot = figure();
xlabel('PER');ylabel('Mean Number of Warnings');hold on;plot(PER_values,meanNumWarnings_camp,'r');legend('Camp_linear');     
hold on;plot(PER_values,meanNumWarnings_knipling,'b');legend('Knipling');
hold on;plot(PER_values,meanNumWarnings_NHTSA_early,'m');plot(PER_values,meanNumWarnings_NHTSA_inter,'g');plot(PER_values,meanNumWarnings_NHTSA_immi,'r');legend('Camp linear', 'Knipling','NHSTA early','NHSTA inter','NHSTA immi');                

Ttitle='Num_Warnings_vs_PER';
path=['/home/nitish/Documents/FCWProject/MatlabCode/' Ttitle '.png'];
saveas(warnings_plot,path);

save FCW_Model_Results1.mat



    %load FCW_Model_Results1.mat
end

%%for Data Communicated over a network, using periodic beaconing
function [v_lv_rcvd,a_lv_rcvd,x_lv_rcvd,PBRATE] = findReceivedLVInfoOverPBNetwork(v_lv,a_lv,x_lv, rate,PER)
%base rate is 10Hz
global CONST_SPEED_MODEL;
numsent = 1;
v_lv_rcvd = v_lv;a_lv_rcvd = a_lv;x_lv_rcvd=x_lv;
if (rate > 10)
    disp('Error, RATE should be less than 10');
    return;
end
interval = 10/rate;lastrxindx=1;datatimestep = (1/10);
for tt = 2:size(v_lv,2)
    rtt = ceil(interval*floor(tt/interval)); % time spaced to achieve the rate
    if(rtt == 0), rtt = tt;end;
    v_lv_rcvd(tt) = v_lv( rtt ); %use sample and hold -> leads to constant speed coasting
    a_lv_rcvd(tt) = a_lv( rtt );
    x_lv_rcvd(tt) = x_lv( rtt );
    lost = (rand(1)<PER) ;%getPERfromNetworkModel();
    if(lost == 0 && tt == rtt) % a new packet received
        lastrxindx = tt;
        numsent = numsent+1;
    end
    if CONST_SPEED_MODEL
        x_lv_rcvd(tt) =  x_lv(lastrxindx) + datatimestep*(tt - lastrxindx)*v_lv(lastrxindx);
        v_lv_rcvd(tt) =  v_lv(lastrxindx) ;  %use sample and hold -> leads to constant speed coasting
        if(tt == lastrxindx)
            a_lv_rcvd(tt) =  a_lv(tt);
        else
            a_lv_rcvd(tt) =  0;
        end
    else
        x_lv_rcvd(tt) =  x_lv(lastrxindx) + datatimestep*(tt - lastrxindx)*v_lv(lastrxindx) + a_lv(lastrxindx)*(datatimestep*(tt - lastrxindx))^2;
        v_lv_rcvd(tt) =  v_lv(lastrxindx) + (tt - lastrxindx)*a_lv(lastrxindx)*datatimestep;  %use sample and hold -> leads to constant speed coasting
        a_lv_rcvd(tt) =  a_lv(lastrxindx) ;
    end
    
end
PBRATE = (datatimestep/0.1)*numsent/(tt*datatimestep);
end

function [r_w] = CAMPLinearWarningAlg(Vlv,Vfv,Alv,Afv) % r is range,
%td is the total delay time including driverresponse time and brake system delay
td=2.05;
G = 9.8;

for tt = 1:size(Vlv,2)
    %predicted SV speed after the delay
    Vfvp(tt) = Vfv(tt) + Afv(tt)*(td); %Asv should be negative
    Vlvp(tt) = Vlv(tt) + Alv(tt)*(td); %Asv should be negative
    if Vfvp(tt) < 0
        Vfvp(tt) = 0;
    end
    if Vlvp(tt) < 0
        Vlvp(tt) = 0;
    end
    
    %required deceleration: ( in ft/s and ft/s^2 )
    %convert m/s to ft/s
    cnv = 3.28084; % m/s to MPH : 2.23694
    dec_req(tt)= -5.308 + 0.685*Alv(tt)*cnv + 2.57*(Vlv(tt)*cnv>0) - 0.086*(Vfvp(tt) - Vlvp(tt))*cnv;
    dec_req(tt) = dec_req(tt) / cnv;    %convert back to m
    
    %the following is in m/s and is correct and similar to the above required decel. (g) = -0.164 + 0.668(lead decel. in g�s) -0.00368(closing speed in MPH) + 0.078(if lead moving)
    %dec_req(tt)= G*( -0.165 + 0.685*(Alv(tt)/G) + 0.080*(Vlv(tt)>0) - 0.00877*(Vfvp(tt) - Vlvp(tt)) );
    %range_lost during td
    r_d(tt) = (Vfv(tt) - Vlv(tt))*td + 0.5*(Afv(tt) - Alv(tt))*(td^2);
    if( Vlv(tt) == 0) % case 1, LV stationary all the way
        BOR = (Vfvp(tt)^2) / (-2*dec_req(tt));
    elseif (Vlv(tt) > 0 && Vlvp(tt) > 0) % Case 2 LV moving , still moving
        dec_lv(tt) = Alv(tt);
        BOR = ((Vfvp(tt) - Vlvp(tt))^2) / (-2*(dec_req(tt) - dec_lv(tt)));
    elseif (Vlv(tt) > 0 && Vlvp(tt) <= 0)
        dec_lv(tt) = Alv(tt);
        BOR = (Vfvp(tt)^2)/(-2*(dec_req(tt))) -  (Vlvp(tt)^2)/(-2*(dec_lv(tt)));
    end
    
    r_w(tt) = BOR + r_d(tt);
    if(r_w(tt) < 0 )
        r_w(tt) = 0;
    end
end
end

function [rw_knipling,lv_state] = Knipling(Vlv,Vfv,Alv,Afv)
td=2.05;
G = 9.8;
cnv = 3.28084;
rw_knipling = zeros(1,size(Vlv,2));
acc_const = 0.6*G;
for tt = 1:size(Vlv,2)
%for tt = 1:20
    %disp(size(Vlv,2));
    
%     Vlvp = cnv*Vlv(tt);
%     Vfvp = cnv*Vfv(tt);
%     Alvp = cnv*Alv(tt);
%     Afvp = cnv*Afv(tt);
    if(Vlv(tt) == 0) % case 1, LV stationary
        %disp('Yes');
        %if(Afv(tt) ~= 0)
            lv_state(tt)=0;
            rw_knipling(tt) = td*Vfv(tt) + (Vfv(tt)^2)/(2*acc_const);
        %else
        %    rw_knipling(tt) = 100;   %50 is just a number to show infinity rw
        %end
    else % Case 2 LV moving
        %disp('No');
        %if(Afv(tt)~=0 && Alv(tt)~=0)
            lv_state(tt)=1;
            rw_knipling(tt) = (Vfv(tt)^2)/(2*acc_const) + (td*Vfv(tt)) + (Vlv(tt)^2)/(2*Alv(tt));
        %else
        %    rw_knipling(tt)=100;
        %end
    end
    
    if(rw_knipling(tt) < 0)
        rw_knipling(tt)=0;
    end
   if(rw_knipling(tt) > 100)
    rw_knipling(tt)=100;
   end    
end
%rw_knipling

%lv_state
end

function [Dmiss_early,Dmiss_inter,Dmiss_immi] = NHTSA(Xlv,Xfv,VehLen,Vlv,Vfv,Alv,Afv)
G = 9.8;
tr_b = 0.5;
tr_nb = 1.6;
Afvmax_early = 0.38*G;
Afvmax_inter = 0.45*G;
Afvmax_immi = 0.55*G;
% early_warning_counter = 0;
% inter_warning_counter = 0;
% immi_warning_counter = 0;
% no_warning_counter = 0;
A_temp_deno = 0.001;
%rw_NHTSA_early = zeros(1,size(Vlv,2));
%rw_NHTSA_inter = zeros(1,size(Vlv,2));
%rw_NHTSA_immi = zeros(1,size(Vlv,2));
Dmiss_early = zeros(1,size(Vlv,2));
Dmiss_inter = zeros(1,size(Vlv,2));
Dmiss_immi = zeros(1,size(Vlv,2));

for tt = 1:size(Vlv,2)
    Vr = Vlv(tt)-Vfv(tt);
    %set tr based on whether the fv is braking or not
    if(Afv(tt) < 0)
        tr = tr_b;
    else
        tr = tr_nb;
    end
    %r = Xlv(tt)-Xfv(tt)-VehLen;    %assuming that r refers to distacne b/w FV and LV
    %r=tr*Vfv(tt);  %assuming that r refers to distance covered during delay
    r = 61; %assuming that r refers to sensor's range
    
    %First step: find time-to-stop of both LV and FV
    if(Alv(tt)<0.001 && Alv(tt)>-0.001)
        tlvs = -Vlv(tt)/A_temp_deno;
    else
        tlvs = -Vlv(tt)/Alv(tt);
    end
    
    if(Vfv(tt)-Afv(tt)*tr >= 0) %if fv doesn't stop before tr is reached
        tfvs_early = tr - ((Vfv(tt)-Afv(tt)*tr)/Afvmax_early);
        tfvs_inter = tr - ((Vfv(tt)-Afv(tt)*tr)/Afvmax_inter);
        tfvs_immi = tr - ((Vfv(tt)-Afv(tt)*tr)/Afvmax_immi);
    else                        %if fv stops before tr is reached
        if(Afv(tt)<0.001 && Afv(tt)>-0.001)
            tfvs_early = -Vfv(tt)/A_temp_deno;
            tfvs_inter = -Vfv(tt)/A_temp_deno;
            tfvs_immi = -Vfv(tt)/A_temp_deno;
        else
            tfvs_early = -Vfv(tt)/Afv(tt);
            tfvs_inter = -Vfv(tt)/Afv(tt);
            tfvs_immi = -Vfv(tt)/Afv(tt);
        end
    end
    
    %Second Step: Find Dmiss
    
    %Calculating Dmiss_early
    if(Vlv(tt)>0 && tlvs<tfvs_early)
        Dmiss_early(tt) = r + ((Afv(tt)-Afvmax_early)*(tr^2))/2 - (Alv(tt)*(tlvs^2))/2 - (Afv(tt)-Afvmax_early)*tr*tfvs_early + Vr*tfvs_early + Alv(tt)*tfvs_early*tlvs - (Afvmax_early*(tfvs_early^2))/2;
    elseif(Vlv(tt)==0 || tfvs_early<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_early-A_temp_deno)) + tr;
            Dmiss_early(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_early)*(tm^2))/2 - (Afv(tt)-Afvmax_early)*tm*tr + ((Afv(tt)-Afvmax_early)*(tr^2))/2;
        else
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_early-Alv(tt))) + tr;
            Dmiss_early(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_early)*(tm^2))/2 - (Afv(tt)-Afvmax_early)*tm*tr + ((Afv(tt)-Afvmax_early)*(tr^2))/2;
        end
    end
    
    %Calculating Dmiss_inter
    if(Vlv(tt)>0 && tlvs<tfvs_inter)
        Dmiss_inter(tt) = r + ((Afv(tt)-Afvmax_inter)*(tr^2))/2 - (Alv(tt)*(tlvs^2))/2 - (Afv(tt)-Afvmax_inter)*tr*tfvs_inter + Vr*tfvs_inter + Alv(tt)*tfvs_inter*tlvs - (Afvmax_inter*(tfvs_inter^2))/2;
    elseif(Vlv(tt)==0 || tfvs_inter<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_inter-A_temp_deno)) + tr;
            Dmiss_inter(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_inter)*(tm^2))/2 - (Afv(tt)-Afvmax_inter)*tm*tr + ((Afv(tt)-Afvmax_inter)*(tr^2))/2;
        else
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_inter-Alv(tt))) + tr;
            Dmiss_inter(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_inter)*(tm^2))/2 - (Afv(tt)-Afvmax_inter)*tm*tr + ((Afv(tt)-Afvmax_inter)*(tr^2))/2;
        end
    end
    
    %Calculating Dmiss_immi
    if(Vlv(tt)>0 && tlvs<tfvs_immi)
        Dmiss_immi(tt) = r + ((Afv(tt)-Afvmax_immi)*(tr^2))/2 - (Alv(tt)*(tlvs^2))/2 - (Afv(tt)-Afvmax_immi)*tr*tfvs_immi + Vr*tfvs_immi + Alv(tt)*tfvs_immi*tlvs - (Afvmax_immi*(tfvs_immi^2))/2;
    elseif(Vlv(tt)==0 || tfvs_immi<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_immi-A_temp_deno)) + tr;
            Dmiss_immi(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_immi)*(tm^2))/2 - (Afv(tt)-Afvmax_immi)*tm*tr + ((Afv(tt)-Afvmax_immi)*(tr^2))/2;
        else
            tm = ((Vr+(Alv(tt)-Afv(tt))*tr)/(Afvmax_immi-Alv(tt))) + tr;
            Dmiss_immi(tt) = r + Vr*tm + ((Alv(tt)-Afvmax_immi)*(tm^2))/2 - (Afv(tt)-Afvmax_immi)*tm*tr + ((Afv(tt)-Afvmax_immi)*(tr^2))/2;
        end
    end
    
    %Third step: Choose the appropriate alert based on comparing d_miss to-
    % the d_threshold
    
%    rw_NHTSA(tt)=Dmiss_immi(tt);
%     Dthresh = r+Vfv(tt)*0.1;

%     Dthresh = r+Vfv(tt)*0.1; %if the threshold is supposed to the same then i can just this out of this func
%     if(Dmiss_immi(tt) < Dthresh)
%         rw_NHTSA(tt) = Dmiss_immi(tt);
%         immi_warning_counter = immi_warning_counter+1;
%         %disp('Early Warning!');
%     elseif(Dmiss_inter(tt) < Dthresh)
%         rw_NHTSA(tt) = Dmiss_inter(tt);
%         inter_warning_counter = inter_warning_counter+1;
%         %disp('Intermediate Warning!!');
%     elseif(Dmiss_early(tt) < Dthresh)
%         early_warning_counter = early_warning_counter+1;
%         %disp('Imminent Warning!!!');
%         rw_NHTSA(tt) = Dmiss_early(tt);
%     else
%         rw_NHTSA(tt) = 0;
%         no_warning_counter = no_warning_counter+1;
%         %disp('No Warning');
%     end

%     if(rw_NHTSA(tt)<0)
%         rw_NHTSA(tt) = 0;
%     end
    if(Dmiss_early(tt)<0)
        Dmiss_early(tt)=0;
    end
    if(Dmiss_inter(tt)<0)
        Dmiss_inter(tt)=0;
    end
    if(Dmiss_immi(tt)<0)
        Dmiss_immi(tt)=0;
    end
    
end
end


function [Dmiss_early,Dmiss_inter,Dmiss_immi] = NHTSA_2(Xlv,Xfv,VehLen,Vlv,Vfv,Alv,Afv)
G = 9.8;
tr_b = 0.5;
tr_nb = 1.6;
Afvmax_early = 0.38*G;
Afvmax_inter = 0.45*G;
Afvmax_immi = 0.55*G;
% early_warning_counter = 0;
% inter_warning_counter = 0;
% immi_warning_counter = 0;
% no_warning_counter = 0;
A_temp_deno = 0.001;
%rw_NHTSA_early = zeros(1,size(Vlv,2));
%rw_NHTSA_inter = zeros(1,size(Vlv,2));
%rw_NHTSA_immi = zeros(1,size(Vlv,2));
Dmiss_early = zeros(1,size(Vlv,2));
Dmiss_inter = zeros(1,size(Vlv,2));
Dmiss_immi = zeros(1,size(Vlv,2));

for tt = 1:size(Vlv,2)
    Vr = Vfv(tt)-Vlv(tt);
    %set tr based on whether the fv is braking or not
    if(Afv(tt) < 0)
        tr = tr_b;
    else
        tr = tr_nb;
    end
    %r = Xlv(tt)-Xfv(tt)-VehLen;    %assuming that r refers to distacne b/w FV and LV
    %r=tr*Vfv(tt);  %assuming that r refers to distance covered during delay
    %r = 61; %assuming that r refers to sensor's range
    r=0.1*Vfv(tt)+Xfv(tt)-Xlv(tt)+VehLen;   %using this from the other paper
    
    %First step: find time-to-stop of both LV and FV
    if(Alv(tt)<0.001 && Alv(tt)>-0.001)
        tlvs = -Vlv(tt)/A_temp_deno;
    else
        tlvs = -Vlv(tt)/Alv(tt);
    end
    
    if(Vfv(tt)-Afv(tt)*tr >= 0) %if fv doesn't stop before tr is reached
        tfvs_early = tr - ((Vfv(tt)+Afv(tt)*tr)/Afvmax_early);
        tfvs_inter = tr - ((Vfv(tt)+Afv(tt)*tr)/Afvmax_inter);
        tfvs_immi = tr - ((Vfv(tt)+Afv(tt)*tr)/Afvmax_immi);
    else                        %if fv stops before tr is reached
        if(Afv(tt)<0.001 && Afv(tt)>-0.001)
            tfvs_early = -Vfv(tt)/A_temp_deno;
            tfvs_inter = -Vfv(tt)/A_temp_deno;
            tfvs_immi = -Vfv(tt)/A_temp_deno;
        else
            tfvs_early = -Vfv(tt)/Afv(tt);
            tfvs_inter = -Vfv(tt)/Afv(tt);
            tfvs_immi = -Vfv(tt)/Afv(tt);
        end
    end
    
    %Second Step: Find Dmiss
    
    %Calculating Dmiss_early

    if(Vlv(tt)>0 && tlvs<tfvs_early)
        Dmiss_early(tt) = r - ((Afv(tt)-Afvmax_early)*(tr^2))/2 + (Alv(tt)*(tlvs^2))/2 + (Afv(tt)-Afvmax_early)*tr*tfvs_early - Vr*tfvs_early - Alv(tt)*tfvs_early*tlvs + (Afvmax_early*(tfvs_early^2))/2;
    elseif(Vlv(tt)==0 || tfvs_early<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr-(Afv(tt)-Alv(tt))*tr)/(Afvmax_early-A_temp_deno)) + tr;
            Dmiss_early(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_early)*(tm^2))/2 + (Afv(tt)-Afvmax_early)*tm*tr - ((Afv(tt)-Afvmax_early)*(tr^2))/2;
        else
            tm = ((Vr-(Afv(tt)-Alv(tt))*tr)/(Afvmax_early-Alv(tt))) + tr;
            Dmiss_early(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_early)*(tm^2))/2 + (Afv(tt)-Afvmax_early)*tm*tr - ((Afv(tt)-Afvmax_early)*(tr^2))/2;
        end
    end
    
    %Calculating Dmiss_inter
    if(Vlv(tt)>0 && tlvs<tfvs_inter)
        Dmiss_inter(tt) = r - ((Afv(tt)-Afvmax_inter)*(tr^2))/2 + (Alv(tt)*(tlvs^2))/2 + (Afv(tt)-Afvmax_inter)*tr*tfvs_inter - Vr*tfvs_inter - Alv(tt)*tfvs_inter*tlvs + (Afvmax_inter*(tfvs_inter^2))/2;
    elseif(Vlv(tt)==0 || tfvs_inter<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr-(Afv(tt)-Alv(tt))*tr)/(Afvmax_inter-A_temp_deno)) + tr;
            Dmiss_inter(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_inter)*(tm^2))/2 + (Afv(tt)-Afvmax_inter)*tm*tr - ((Afv(tt)-Afvmax_inter)*(tr^2))/2;
        else
            tm = ((Vr+(Afv(tt)-Alv(tt))*tr)/(Afvmax_inter-Alv(tt))) + tr;
            Dmiss_inter(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_inter)*(tm^2))/2 + (Afv(tt)-Afvmax_inter)*tm*tr - ((Afv(tt)-Afvmax_inter)*(tr^2))/2;
        end
    end
    
    %Calculating Dmiss_immi
    if(Vlv(tt)>0 && tlvs<tfvs_immi)
        Dmiss_immi(tt) = r - ((Afv(tt)-Afvmax_immi)*(tr^2))/2 + (Alv(tt)*(tlvs^2))/2 + (Afv(tt)-Afvmax_immi)*tr*tfvs_immi - Vr*tfvs_immi - Alv(tt)*tfvs_immi*tlvs + (Afvmax_immi*(tfvs_immi^2))/2;
    elseif(Vlv(tt)==0 || tfvs_immi<tlvs)
        if(Alv(tt)<0.001 && Alv(tt)>-0.001)
            tm = ((Vr-(Afv(tt)-Alv(tt))*tr)/(Afvmax_immi-A_temp_deno)) + tr;
            Dmiss_immi(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_immi)*(tm^2))/2 + (Afv(tt)-Afvmax_immi)*tm*tr - ((Afv(tt)-Afvmax_immi)*(tr^2))/2;
        else
            tm = ((Vr-(Afv(tt)-Alv(tt))*tr)/(Afvmax_immi-Alv(tt))) + tr;
            Dmiss_immi(tt) = r - Vr*tm - ((Alv(tt)-Afvmax_immi)*(tm^2))/2 + (Afv(tt)-Afvmax_immi)*tm*tr - ((Afv(tt)-Afvmax_immi)*(tr^2))/2;
        end
    end
    
    %Third step: Choose the appropriate alert based on comparing d_miss to-
    % the d_threshold
    
%    rw_NHTSA(tt)=Dmiss_immi(tt);
%     Dthresh = r+Vfv(tt)*0.1;

%     Dthresh = r+Vfv(tt)*0.1; %if the threshold is supposed to the same then i can just this out of this func
%     if(Dmiss_immi(tt) < Dthresh)
%         rw_NHTSA(tt) = Dmiss_immi(tt);
%         immi_warning_counter = immi_warning_counter+1;
%         %disp('Early Warning!');
%     elseif(Dmiss_inter(tt) < Dthresh)
%         rw_NHTSA(tt) = Dmiss_inter(tt);
%         inter_warning_counter = inter_warning_counter+1;
%         %disp('Intermediate Warning!!');
%     elseif(Dmiss_early(tt) < Dthresh)
%         early_warning_counter = early_warning_counter+1;
%         %disp('Imminent Warning!!!');
%         rw_NHTSA(tt) = Dmiss_early(tt);
%     else
%         rw_NHTSA(tt) = 0;
%         no_warning_counter = no_warning_counter+1;
%         %disp('No Warning');
%     end

%     if(rw_NHTSA(tt)<0)
%         rw_NHTSA(tt) = 0;
%     end
    if(Dmiss_early(tt)<0)
        Dmiss_early(tt)=0;
    end
    if(Dmiss_inter(tt)<0)
        Dmiss_inter(tt)=0;
    end
    if(Dmiss_immi(tt)<0)
        Dmiss_immi(tt)=0;
    end
    
    if(Dmiss_early(tt)>100)
        Dmiss_early(tt)=100;
    end
    if(Dmiss_inter(tt)>100)
        Dmiss_inter(tt)=100;
    end
    if(Dmiss_immi(tt)>100)
        Dmiss_immi(tt)=100;
    end

    
end
end


function [Acc,Prc,TrPos,TrNeg,FlsPos,FlsNeg,gMeanTPP] = analyzeaccuracy(PredictedWrnRange, GroundTruthRw, est_separation_dist, GroundTruthSepDist)
% a,b,c,d from [Lee, Peng] UMichigan paper confusion matrix
% Confusion matrix.
%                             Actual data
%     pred                 Negative   Positive
% Negative (safe)               a       c
% Positive (threatening)        b       d
% Here we limit the search to areas where the distance between cars is in
% terms of time to crash (TODO) is less than 20m

importantRange = 20;
distancetowarning = GroundTruthSepDist - GroundTruthRw;

Acc = 0;Prc = 0;TrPos=0; TrNeg = 0; FlsPos= 0; FlsNeg = 0;

GTThreat = (GroundTruthRw>=GroundTruthSepDist);
PredThreat = (PredictedWrnRange>=est_separation_dist);
GTSafe = (GroundTruthRw<GroundTruthSepDist);
PredSafe = (PredictedWrnRange < est_separation_dist);
a = sum((GTSafe == PredSafe) & (GTSafe == 1) & (distancetowarning < importantRange));
d = sum(GTThreat == PredThreat & (GTThreat == 1) & (distancetowarning < importantRange));
b = sum((GTSafe == PredThreat) & (GTSafe == 1) & (distancetowarning < importantRange));
c = sum((GTThreat == PredSafe) & (GTThreat == 1) & (distancetowarning < importantRange)) ;

if( a+d > 0)
    Acc = (a+d)/(a+b+c+d);
end

if ( b+d > 0)
    Prc = d/(b+d);
end

if( c+d>0)
    TrPos = d/(c+d);
    FlsNeg = c/(c+d);
end

if( a+b>0 )
    TrNeg = a/(a+b);
    FlsPos= b/(a+b);
end

gMeanTPP = Acc;%(TrPos*Prc)^0.5;
% gMeanTPP = (FlsNeg*FlsPos)^0.5;

if(isnan(gMeanTPP) )
    msg = 'NaN found'
end


end

function [Acc,Prc,TrPos,TrNeg,FlsPos,FlsNeg,gMeanTPP] = analyzeaccuracy_NHTSA(threat_vec, safe_vec, GroundTruthRw, est_separation_dist, GroundTruthSepDist)
% a,b,c,d from [Lee, Peng] UMichigan paper confusion matrix
% Confusion matrix.
%                             Actual data
%     pred                 Negative   Positive
% Negative (safe)               a       c
% Positive (threatening)        b       d
% Here we limit the search to areas where the distance between cars is in
% terms of time to crash (TODO) is less than 20m

importantRange = 20;
distancetowarning = GroundTruthSepDist - GroundTruthRw;

Acc = 0;Prc = 0;TrPos=0; TrNeg = 0; FlsPos= 0; FlsNeg = 0;

GTThreat = (GroundTruthRw>=GroundTruthSepDist);
PredThreat = threat_vec;
GTSafe = (GroundTruthRw<GroundTruthSepDist);
PredSafe = safe_vec;

a = sum((GTSafe == PredSafe) & (GTSafe == 1) & (distancetowarning < importantRange));
d = sum(GTThreat == PredThreat & (GTThreat == 1) & (distancetowarning < importantRange));
b = sum((GTSafe == PredThreat) & (GTSafe == 1) & (distancetowarning < importantRange));
c = sum((GTThreat == PredSafe) & (GTThreat == 1) & (distancetowarning < importantRange));

if( a+d > 0)
    Acc = (a+d)/(a+b+c+d);
end

if ( b+d > 0)
    Prc = d/(b+d);
end

if( c+d>0)
    TrPos = d/(c+d);
    FlsNeg = c/(c+d);
end

if( a+b>0 )
    TrNeg = a/(a+b);
    FlsPos= b/(a+b);
end

gMeanTPP = Acc;%(TrPos*Prc)^0.5;
% gMeanTPP = (FlsNeg*FlsPos)^0.5;

if(isnan(gMeanTPP) )
    msg = 'NaN found'
end


end

function [AminusNormal] = aminusf(v)
if (v < 6.1)
    AminusNormal = -8.8;
elseif ( v<12.2)
    AminusNormal = -5.2;
elseif ( v<18.3)
    AminusNormal = -4.4;
elseif ( v<24.4)
    AminusNormal = -2.9;
else
    AminusNormal = -2;
end

end

function [AplusNormal] = aplusf(v)
if (v < 6.1)
    AplusNormal = 7.8;
elseif ( v<12.2)
    AplusNormal = 6.7;
else
    AplusNormal = 4.8;
end

end

function filtered_acc = filterAcc(accel, TAPS)
NN = size(accel,1);
filtered_acc = accel;
timestep = 0.1;
Taprange = floor(TAPS/2);
for t=Taprange+1:NN-Taprange
    filtered_acc(t) = sum(accel(t-Taprange: t+Taprange))/(2*Taprange+1);
    %        est_x(t) = est_x(t-1) + v_lv_est(t-1)*timestep ; % constant speed model
end
end





% BSMs---------------------------------------------------------------------
% Connected V2V safety applications are built around the
% SAE J2735 BSM, which has two parts
% ? BSM Part 1:
% ? Contains the core data elements (vehicle size, position, speed, heading
% acceleration, brake system status)
% ? Transmitted approximately 10x per second
% ? BSM Part 2:
% ? Added to part 1 depending upon events (e.g., ABS activated)
% ? Contains a variable set of data elements drawn from many optional
% data elements (availability by vehicle model varies)
% ? Transmitted less frequently
% ? No on-vehicle BSM storage of BSM data
% ? The BSM is transmitted over DSRC (range ~1,000 meters)

% probTx = 1 - exp(-alpha*e(t-1));
% recvd = (rand<probTx);
%
% estimation every T second
% if recvd == 1
%     txCnt = txCnt+1;

%--------------------------------------------------------------------------
% automotive Acceleration (g)
% event 	typical car 	sports car 	F-1 race car 	large truck
% starting 	0.3�0.5 	0.5�0.9 	1.7 	< 0.2
% braking 	0.8�1.0 	1.0�1.3 	2 	~ 0.6
% cornering 	0.7�0.9 	0.9�1.0 	3 	??