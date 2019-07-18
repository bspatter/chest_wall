function varargout=attenuation_comparisons_script()
clearvars; close all; clc
total_timer=tic;

% load('usData.mat','usData')
%% Ultrasound Attenuation Calculations
% Compares two recorded acoustic signals of a pulse tranmitted through different media (e.g., water and tissue) to determine acoustic attenuation coefficient.

% Import all of the ultrasound data from a table
flags.force_resample=true;
usData= import_us_data(flags.force_resample);

% Tissue parameters
rho_saline09 = 1004.6; % kg/m^3
c_saline09 = 1534; % From Duck, 1990 properties of tissue, Table 4.9. >> interp2([30,40],[0.5;1.0],[1514.4, 1533.9; 1519.7, 1538.9],38,0.9)

rho_tissue = 1000;
c_tissue = 1540; % m/s

rho1 = rho_saline09;
c1 = c_saline09;
z1 = rho1*c1;

rho2 = rho_tissue;
c2 = c_tissue;
z2 = rho2*c2;


filter=struct('Power',unique(usData.Power),'SampleNumber',unique(usData.SampleNumber),'usProbe',unique(usData.usProbe));
filter.Power=-0:-2:-6;
filter.SampleNumber=[35702, 35704, 35710, 35711, 35712, 35714, 35715, 35719, 35724, 35726, 35727, 35737, 35739];
filter.usProbe=categorical({'3s','5s'});
filter.ReferenceIndex=unique(usData.ReferenceIndex(usData.ReferenceIndex>=0));

filter_measurements;



%% compare two tek files to get the attenuation as a function of frequency
tissue_SignalNums = usData.Index(usData.Medium=='rib'); % look at all tissue cases
source_frequencies = unique(usData.frequency_MHz);

% preallocate variable before the main loop
attenuation_timer=tic;
fprintf('Calculating attenuation coefficients...')
fprintf('\n0%% -- ')

[alpha_f,FrequencyBand_]=deal(cell(size(tissue_SignalNums)));
[alpha_fc,fc_,alpha_fmax,alpha_fcs,alpha_dbcm]=deal(zeros(size(tissue_SignalNums)));
for k = 1:numel(tissue_SignalNums)
    tissue_SignalNum=tissue_SignalNums(k);
    tissue_SignalNumi = find(usData.Index==tissue_SignalNum);    
            
    thickness = usData.TissueThickness(tissue_SignalNumi);
    ref_SignalNum = abs(usData.ReferenceIndex(tissue_SignalNumi)); % abs, because - indicates data we will later ignore
    ref_SignalNumi = find(usData.Index==ref_SignalNum);
    
%   Get the time, voltage, and pressure from the recorded tissue signal and its reference
    time1 = usData.Time(ref_SignalNumi,:);
    voltage1 = usData.Voltage(ref_SignalNumi,:);
    pressure1= usData.Pressure(ref_SignalNumi,:);
    
    time2 = usData.Time(tissue_SignalNumi,:);
    voltage2 = usData.Voltage(tissue_SignalNumi,:);
    pressure2= usData.Pressure(tissue_SignalNumi,:);
        
    % Align the peak min pressures in time and start the minimum time at zero
    [~,pmini1]=min(pressure1);
    time1=time1-time1(pmini1);
    [~,pmini2]=min(pressure2);
    time2=time2-time2(pmini2);
    
    tmin = min([time1(:);time2(:)]);
    time1 = time1-tmin;
    time2 = time2-tmin;
    
%     [time1, voltage1, pressure1,time2,voltage2,pressure2] = get_clean_signals(tissue_SignalNumi,ref_SignalNumi);
    
    % Get the FFT of the pressure signals.
    [~,FrequencySpectrum1,pressureFFT1,fc1,pulseBandwidth1]  = find_pulse_frequency(time1,pressure1);
    [~,FrequencySpectrum2,pressureFFT2,fc2,~]  = find_pulse_frequency(time2,pressure2);
    pressureFFT1=reshape(pressureFFT1,[],1);
    pressureFFT2=reshape(pressureFFT2,[],1);
    f1 = pulseBandwidth1*-0.5+fc1;
    f2 = pulseBandwidth1*0.5+fc1;
    pulseBandwidth_index = FrequencySpectrum1>=f1 & FrequencySpectrum1<=f2;
       
    % Calculate alpha
    % calculate_alpha = @(z1,z2,d,A2,A1) 20/d*log10(4*z1*z2/(z1+z2)^2)-20/d*log10(A2./A1);
    depth = round(thickness,3)*100; % Take to nearest mm, convert to cm.
    alpha_f = 20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(abs(pressureFFT2)./abs(pressureFFT1)); % alpha comparing signals 1 and 4

    % Extract frequency and alpha for plotting
     
    fband= FrequencySpectrum1(pulseBandwidth_index)/1e6;
    alpha_dbcm = alpha_f(pulseBandwidth_index);
    alpha_fband = alpha_f(pulseBandwidth_index)./fband;
    FrequencyBand_{k}=fband;
    alpha_{k}=alpha_fband;
    fc_(k)=fc1/1e6;
    alpha_fc(k)=interp1(fband,alpha_fband,fc1/1e6);
    alpha_fmax(k)=(20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(max(abs(pressureFFT2))./max(abs(pressureFFT1))))/((fc1+fc2)/2/1e6);
    alpha_fcs(k)=(20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(interp1(FrequencySpectrum2,abs(pressureFFT2),fc2))./interp1(FrequencySpectrum1,abs(pressureFFT1),fc1))/((fc1+fc2)/2/1e6);
    
    alpha_dbcm_fc(k)=interp1(fband,alpha_dbcm,fc1/1e6);
    
    % Add the attenuation data back into the US Data structure    
    usData.alpha_fc(tissue_SignalNumi)=alpha_fc(k);
    usData.alpha_fmax(tissue_SignalNumi)=alpha_fmax(k);
    usData.alpha_fcs(tissue_SignalNumi)=alpha_fcs(k);
    usData.alpha_dbcm_fc(tissue_SignalNumi)=alpha_dbcm_fc(k);


    % print progress
    myprogress = round(k/numel(tissue_SignalNums)*100);
    if any(k== round(numel(tissue_SignalNums)*(0:0.1:1)))
        fprintf('%03.0d%% -- ',myprogress)
    end    
end
fprintf('\n')
toc(attenuation_timer)


% Get just the tissue cases (i.e., those w/ relevant attenuation data)
tisData=usData(usData.Medium=='rib' & usData.ReferenceIndex>0,:);

if false
figs.line_dBcmMHz_measurement=figure;
for np=1:length(alpha_)
    plot3(FrequencyBand_{np},tissue_SignalNums(np)*ones(size(alpha_{np})),alpha_{np})
    hold on
end
hold off
% ylim([1,2.5])
axis tight
xlabel('Frequency (MHz)')
ylabel('Sample Index')
zlabel('Attenuation Coefficient \alpha (dB / cm / MHz)' )
grid on; box on;
% view(0,90) % Freq vs sample number view
view(0,0) % Freq vs attenuation view
end

if false
figs.line_dBcm_measurement=figure;
for np=1:length(alpha_)
    plot3(FrequencyBand_{np},tissue_SignalNums(np)*ones(size(alpha_{np})),alpha_{np}.*FrequencyBand_{np})
    hold on
end
hold off
% ylim([1,2.5])
axis tight
xlabel('Frequency (MHz)')
ylabel('Sample Index')
zlabel('Attenuation Coefficient \alpha (dB / cm)' )
grid on; box on;
% view(0,90) % Freq vs sample number view
view(0,0) % Freq vs attenuation view
end

plot_SignalNums = usData.SampleNumber(usData.Medium=='rib'); % look at all tissue cases
plot_SignalPowers = usData.Power(usData.Medium=='rib'); % look at all tissue cases



figs.scat_dBcm_donor=figure('Name','Attenuation (dB/cm)');
gs=gscatter(fc_,alpha_fc.*fc_,plot_SignalNums);
axis tight
xlabel('Frequency (MHz)')
ylabel('Attenuation Coefficient \alpha (dB / cm)' )
[gs.MarkerSize]=deal(25);
grid on; box on; hold on
spiffyp(gcf)
xlim([0, max(fc_)+1])
ylim([0, max(alpha_fc.*fc_)+1])

figs.scat_dBcmMHz_donor=figure;
gs=gscatter(fc_,alpha_fc,plot_SignalNums);
set(findobj(figs.scat_dBcmMHz_donor,'type','legend'),'visible','off');
axis tight
xlabel('Frequency (MHz)')
ylabel('Attenuation Coefficient \alpha (dB / cm / MHz)' )
[gs.MarkerSize]=deal(25);
grid on; box on; hold on
spiffyp(gcf)
xlim([0, max(fc_)+1])
ylim([0, max(alpha_fc)+1])




% regression line for each patient
samplenums = unique(plot_SignalNums);
myfitobjs=cell(length(samplenums),1);
for ns = 1:length(samplenums)
    samplenum = samplenums(ns);
    if samplenum<35702;         continue;     end % skip early data, before I had a higher frequency probe   
    ii=(plot_SignalNums==samplenum);% & usData.Power(tissue_SignalNums)==-6;
    
    xx=fc_(ii);
    yy=alpha_fc(ii);
    
    % remove high frequency set (anything above 6 MHz)
    yy=yy(xx<6);
    xx=xx(xx<6);
    
    p=polyfit(xx,yy,1);
    xp1=linspace(min(xx),max(xx),100);
    yp1=polyval(p,xp1);
    xp2=linspace(0,min(xx));
    yp2=polyval(p,xp2);
    
    plot(xp1,yp1,'Color',gs(ns).Color);
    plot(xp2,yp2,'linestyle','--','Color',gs(ns).Color)
    text(xp1(end)*1.02,yp1(end),sprintf('%.0f',samplenum))
    
    
    xx=fc_(ii);
    yy=alpha_dbcm_fc(ii);
    
    % remove high frequency set (anything above 6 MHz)
    yy=yy(xx<6);
    xx=xx(xx<6);
%     display(samplenum)
    myfitobjs{ns}=fit(xx(:),yy(:),'power2');    
end
myfitobjs(cellfun(@isempty,myfitobjs))=[];
xline(min(xx))
hold off
xlabel('Frequency (MHz)')
ylabel('Attenuation Coefficient \alpha (dB / cm / MHz)' )
title('Best fit \alpha(f_c) lines per donor')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
xlim([1,max(xx)+1])


figs.scat_dBcm_power=figure('Name','Attenuation (dB/cm)');
gs=gscatter(fc_,alpha_fc.*fc_,plot_SignalPowers);
axis tight
xlabel('Frequency (MHz)')
ylabel('Attenuation Coefficient \alpha (dB / cm)' )
[gs.MarkerSize]=deal(25);
grid on; box on; hold on
spiffyp(gcf)
xlim([0, max(fc_)+1])
ylim([0, max(alpha_fc.*fc_)+1])

% Analyze data by source frequency / power
source_frequencies = source_frequencies(source_frequencies<6);
donors = unique(tisData.SampleNumber); donors = donors(donors>35691);
[alpha_fc_donor_mean,alpha_fc_donor_std]=deal(zeros(length(source_frequencies),length(donors)));
[alpha_fc_mean,alpha_fc_std,alpha_dbcm_fc_mean,alpha_dbcm_fc_std]=deal(zeros(1,numel(source_frequencies)));
for ns = 1:length(source_frequencies)
    fii = tisData.frequency_MHz==source_frequencies(ns); % All data w/ the current source frequency
    alpha_fc_mean(ns) = mean(tisData.alpha_fc(fii));
    alpha_fc_std(ns) = std(tisData.alpha_fc(fii));
    
    alpha_dbcm_fc_mean(ns) = mean(tisData.alpha_dbcm_fc(fii));
    alpha_dbcm_fc_std(ns) = std(tisData.alpha_dbcm_fc(fii));
    
    for nd = 1:length(donors)
        dii = tisData.SampleNumber==donors(nd); % All data w/ the current source frequency
        alpha_fc_donor_mean(ns,nd) = mean(tisData.alpha_fc(fii & dii));
        alpha_fc_donor_std(ns,nd) = std(tisData.alpha_fc(fii & dii));        
    end    
end

figs.box_dBcmMHz=figure;
boxplot(tisData.alpha_fc,tisData.frequency_MHz);

if false
% Source frequency vs alpha dB/cm/MHz
figs.line_dBcmMHz_MeanStd=figure; axes; hold on;
plot(source_frequencies, alpha_fc_mean,'-o')
plot(source_frequencies, alpha_fc_mean-alpha_fc_std,'r-o')
plot(source_frequencies, alpha_fc_mean+alpha_fc_std,'r-o')
xlabel('Frequency setting (MHz)')
ylabel('Attenuation coefficient (dB/cm/MHz)')
xlim([0, max(source_frequencies)+1])
ylim([0, alpha_fc_mean(end)+alpha_fc_std(end)+1])
end

figs.line_dBcm_MeanStd=figure; axes; hold on;
plot(source_frequencies, alpha_dbcm_fc_mean,'-o')
plot(source_frequencies, alpha_dbcm_fc_mean-alpha_dbcm_fc_std,'r-o')
plot(source_frequencies, alpha_dbcm_fc_mean+alpha_dbcm_fc_std,'r-o')
xlabel('Frequency setting (MHz)')
ylabel('Attenuation coefficient (dB/cm)')



figs.line_dBcmMHz_donor=figure; axes; hold on; %#ok<*STRNU>
plot(source_frequencies, alpha_fc_donor_mean)
xlabel('Frequency setting (MHz)')
ylabel('Attenuation coefficient (dB/cm/MHz)')





toc(total_timer)
disp('done');













function filter_measurements
    usData=usData(ismember(usData.SampleNumber, filter.SampleNumber),:);
    usData=usData(ismember(usData.Power,        filter.Power),:);
    usData=usData(ismember(usData.usProbe,      filter.usProbe),:);
    usData=usData(ismember(usData.ReferenceIndex, filter.ReferenceIndex),:);
end

function [time1, voltage1, pressure1,time2,voltage2,pressure2] = get_clean_signals(tissue_SignalNumi,ref_SignalNumi) %#ok<DEFNU>
        % Get the time, voltage, and pressure from the recorded tissue signal and its reference
    time1 = usData.Time(ref_SignalNumi,:);
    time1 = time1-time1(1);
    voltage1 = usData.Voltage(ref_SignalNumi,:);
    voltage1=voltage1-mean(voltage1);
    [time1,voltage1]=zpad(time1,voltage1,1e5);
    pressure1 = calibrate_V2Pa(time1,voltage1,'HGL0200_calibration.csv');

    time2 = usData.Time(tissue_SignalNumi,:);
    time2 = time2-time2(1);
    voltage2 = usData.Voltage(tissue_SignalNumi,:);
    voltage2=voltage2-mean(voltage2);
    [time2,voltage2]=zpad(time2,voltage2,1e5);
    pressure2 = calibrate_V2Pa(time2,voltage2,'HGL0200_calibration.csv');

    % Align the peak min pressures in time and start the minimum time at zero
    [~,pmini1]=min(pressure1);
    time1=time1-time1(pmini1);
    [~,pmini2]=min(pressure2);
    time2=time2-time2(pmini2);

    tmin = min([time1(:);time2(:)]);
    time1 = time1-tmin;
    time2 = time2-tmin;
end

% % 
% % function get_attenuation_coeffs
% % attenuation_timer=tic;
% % fprintf('Calculating attenuation coefficients...')
% % fprintf('\n0%% -- ')
% % 
% % % preallocate variable before the main loop
% % [alpha_f,FrequencyBand_]=deal(cell(size(tissue_SignalNums)));
% % [alpha_fc,fc_,alpha_fmax,alpha_fcs,alpha_dbcm]=deal(zeros(size(tissue_SignalNums)));
% % 
% % 
% % 
% % for k = 1:numel(tissue_SignalNums)
% %     tissue_SignalNum=tissue_SignalNums(k);
% %     tissue_SignalNumi = find(usData.Index==tissue_SignalNum);    
% %             
% %     thickness = usData.TissueThickness(tissue_SignalNumi);
% %     ref_SignalNum = usData.ReferenceIndex(tissue_SignalNumi);
% %     ref_SignalNumi = find(usData.Index==ref_SignalNum);
% %     
% % %   Get the time, voltage, and pressure from the recorded tissue signal and its reference
% %     time1 = usData.Time(ref_SignalNumi,:);
% %     voltage1 = usData.Voltage(ref_SignalNumi,:);
% %     pressure1= usData.Pressure(ref_SignalNumi,:);
% %     
% %     time2 = usData.Time(tissue_SignalNumi,:);
% %     voltage2 = usData.Voltage(tissue_SignalNumi,:);
% %     pressure2= usData.Pressure(tissue_SignalNumi,:);
% %         
% %     % Align the peak min pressures in time and start the minimum time at zero
% %     [~,pmini1]=min(pressure1);
% %     time1=time1-time1(pmini1);
% %     [~,pmini2]=min(pressure2);
% %     time2=time2-time2(pmini2);
% %     
% %     tmin = min([time1(:);time2(:)]);
% %     time1 = time1-tmin;
% %     time2 = time2-tmin;   
% %     
% %     % Get the FFT of the pressure signals.
% %     [~,FrequencySpectrum1,pressureFFT1,fc1,~]  = find_pulse_frequency(time1,pressure1);
% %     [~,FrequencySpectrum2,pressureFFT2,fc2,~]  = find_pulse_frequency(time2,pressure2);
% %     pressureFFT1=reshape(pressureFFT1,[],1);
% %     pressureFFT2=reshape(pressureFFT2,[],1);
% %     f1 = pulseBandwidth1*-0.5+fc1;
% %     f2 = pulseBandwidth1*0.5+fc1;
% %     pulseBandwidth_index = FrequencySpectrum1>=f1 & FrequencySpectrum1<=f2;
% %        
% %     % Calculate alpha
% %     % calculate_alpha = @(z1,z2,d,A2,A1) 20/d*log10(4*z1*z2/(z1+z2)^2)-20/d*log10(A2./A1);
% %     depth = round(thickness,3)*100; % Take to nearest mm, convert to cm.
% %     alpha_f = 20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(abs(pressureFFT2)./abs(pressureFFT1)); % alpha comparing signals 1 and 4
% % 
% %     % Extract frequency and alpha for plotting
% %     fband= FrequencySpectrum1(pulseBandwidth_index)/1e6;
% %     alpha_dbcm = alpha_f(pulseBandwidth_index);
% %     alpha_fband = alpha_f(pulseBandwidth_index)./fband;
% % %     FrequencyBand_{k}=fband;
% % %     alpha_{k}=alpha_fband;
% % %     fc_(k)=fc1/1e6;
% % %     alpha_fc(k)=interp1(fband,alpha_fband,fc1/1e6);
% %     usData(tissue_SignalNumi).RefFrequencyBand=fband;
% %     usData(tissue_SignalNumi).alpha_RefFrequencyBand=alpha_fband;
% %     usData(tissue_SignalNumi).ref_fc=fc1/1e6;
% %     alpha_fmax(k)=(20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(max(abs(pressureFFT2))./max(abs(pressureFFT1))))/((fc1+fc2)/2/1e6);
% %     alpha_fcs(k)=(20/depth*log10(4*z1*z2/(z1+z2)^2)-20/depth*log10(interp1(FrequencySpectrum2,abs(pressureFFT2),fc2))./interp1(FrequencySpectrum1,abs(pressureFFT1),fc1))/((fc1+fc2)/2/1e6);
% %     
% %     alpha_dbcm_fc(k)=interp1(fband,alpha_dbcm,fc1/1e6);
% %     
% %     % Add the attenuation data back into the US Data structure    
% %     usData.alpha_fc(tissue_SignalNumi)=alpha_fc(k);
% %     usData.alpha_fmax(tissue_SignalNumi)=alpha_fmax(k);
% %     usData.alpha_fcs(tissue_SignalNumi)=alpha_fcs(k);
% %     usData.alpha_dbcm_fc(tissue_SignalNumi)=alpha_dbcm_fc(k);
% % 
% %     % print progress
% %     myprogress = round(k/numel(tissue_SignalNums)*100);
% %     if any(k== round(numel(tissue_SignalNums)*(0:0.1:1)))
% %         fprintf('%03.0d%% -- ',myprogress)
% %     end    
% % end
% % fprintf('\n')
% % toc(attenuation_timer)
% % 
% % 
% % 
% % 
% % end

end



% ToDo
%%
% * Calculate Mean and STD of attenuation as a function of frequency
% Calculate a regression line for each patient
%
% Done
% * Add reflection (Note: very small part of alpha, so density doesn't matter
% much).
% * make choice of depth/reference file automatic.
% * Loop over various data points and combine plots of alpha in bandwith













%%% OLD PIECES
% Get the FFT of the pressure signal.
% time_adjusted1 = time1-time1(1); % start time from 0;
% TimeTotal1=range(time_adjusted1); % Total time
% SamplingFrequency1 = 1/mean(diff(time_adjusted1)); % sampling frequency
% N1 = round(TimeTotal1*SamplingFrequency1); %number of data points
% N1 = round(N1/2)*2; %force it to be even
% TimeTotal1 = N1/SamplingFrequency1; %new total time domain length, if N changed
% TimeFFT1 = (0:(N1-1))/SamplingFrequency1; % time vector (max(t) is 1 bin short of T)
% FrequencySpectrum1 = (0:(N1/2))'/N1*SamplingFrequency1; %frequency vector (min(f)=0, the DC component. max(f) = fs/2 exactly)
% % Match pressure to output variables
% if length(TimeFFT1)==length(pressure1)
%     y = pressure1;%function(t); %hypothetical time domain vector, length(y)=N
% else
%     y = interp1(time_adjusted1,pressure1, TimeFFT1);%function(t); %hypothetical time domain vector, length(y)=N
% end
% % lpad=2*length(time); % Do a bit of zero padding.
% pressure1FFT= bfft(y);


% Get the FFT of the pressure signal.
% time_adjusted2 = time2-time2(1); % start time from 0;
% TimeTotal2=range(time_adjusted2); % Total time
% SamplingFrequency2 = 1/mean(diff(time_adjusted2)); % sampling frequency
% N2 = round(TimeTotal2*SamplingFrequency2); %number of data points
% N2 = round(N2/2)*2; %force it to be even
% TimeTotal2 = N2/SamplingFrequency2; %new total time domain length, if N changed
% TimeFFT2 = (0:(N2-1))/SamplingFrequency2; % time vector (max(t) is 1 bin short of T)
% FrequencySpectrum2 = (0:(N2/2))'/N2*SamplingFrequency2; %frequency vector (min(f)=0, the DC component. max(f) = fs/2 exactly)
% % Match pressure to output variables
% if length(TimeFFT2)==length(pressure2)
%     y = pressure2;%function(t); %hypothetical time domain vector, length(y)=N
% else
%     y = interp1(time_adjusted2,pressure2, TimeFFT2);%function(t); %hypothetical time domain vector, length(y)=N
% end
% % lpad=2*length(time); % Do a bit of zero padding.
% pressure2FFT= bfft(y);