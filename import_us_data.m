function A = import_us_data(varargin)

flags.force_resample=false;
if numel(varargin)>0
    if varargin{1}==true
        flags.force_resample=true;
    end
end


% import base table
A=readtable('measurement_table.xlsx');
B=readtable('measurement_table_old.xlsx');

if isequal(A,B) && ~flags.force_resample
    load('usData.mat','usData');
    A=usData;
else
    
    % modify table for processing
    A.usProbe=categorical(A.usProbe);
    A.Medium=categorical(A.Medium);
    A.Hydrophone=categorical(A.Hydrophone);
    A.date=datetime(A.date);
    A=addvars(A, (1:size(A,1))','Before','SampleNumber','NewVariableNames','Index');
    A=addvars(A, zeros(size(A,1),1),'After','frequency_MHz','NewVariableNames','alpha_fc');
    
    %% Read the tek files to get time, voltage. Calcualte pressure and reference indices
    get_time_voltage;
    get_reference_indices;
    get_pressure;
    
    
    %% Save the output
    usData=A;
    save('usData.mat','usData');
    
    %Copy the file to look for differences next time
    copyfile('measurement_table.xlsx','measurement_table_old.xlsx')
end


% Zero pad the time and voltage matrices to prepare for FFT
    function [Time,Voltage]=zeropad_tv(Time,Voltage)
        %         Npad=1e5;
        Npad=2^nextpow2(size(Time,2));
        Time_pad = Time(:,end)+mean(diff(Time,1,2),2)*(1:(Npad-size(Time,2)));
        Time = [Time, Time_pad];
        
        Voltage = [Voltage, zeros(size(Voltage,1), Npad-size(Voltage,2) )];
    end

% Calculate the pressure using the time, voltage, and hydrophone senstivity curve
    function get_time_voltage
        tv_timer=tic;
        fprintf('Getting time and voltage data from .tek files...\n')
        %preallocate
        [Time,Voltage] = deal(zeros(size(A,1),1e4));
        
        for n=(A.Index)'
            fname=sprintf("E:/brandon/research/chest_wall_attenuation/wavedata/human_ribs/s%.0f-%s/tek%04.0fCH1.csv",A.SampleNumber(n),A.usProbe(n),A.tekNumber(n));
            [t,v]=read_input_data(fname);
            SamplingFrequency=1/diff(t(1:2));
%             v = lowpass(v,100e6,SamplingFrequency); % Cleans the signals, but slows things down a lot and does not appear to make much difference
            Time(n,:)=reshape(t,1,[]);
            Voltage(n,:)=reshape(v,1,[]);
        end
        A.Time = Time;
        
        A.Voltage = Voltage;
        toc(tv_timer);
        fprintf('\n')
    end

% % Calculate the pressure using the time, voltage, and hydrophone senstivity curve
%     function get_pressure
%         pressure_timer=tic;
%         fprintf('Converting voltages to pressure...\n')
%         % Adjust the time to start at zero subtract the mean voltage (must be done after get_reference_indices)
%         A.Time = A.Time-A.Time(:,1);
%         A.Voltage = A.Voltage-mean(A.Voltage,2);
%         
%         [A.Time,A.Voltage]=zeropad_tv(A.Time,A.Voltage);
%         
%         p1=zeros(size(A.Time));
%         fprintf('0%% -- ')
%         calibration_data=readtable('HGL0200_calibration.csv');
% %         calibration_data='HGL0200_calibration.csv';
%         for n = 1:size(A.Time,1)
%             p1(n,:) = reshape(calibrate_V2Pa(A.Time(n,:),A.Voltage(n,:),calibration_data),1,[]);
%             
%             % print progress
%             myprogress = round(n/size(A.Time,1)*100);
%             if any(n== round(size(A.Time,1)*(0:0.1:1)))
%                 fprintf('%03.0d%% -- ',myprogress)
%             end
%         end
%         A=addvars(A, p1,'After','Time','NewVariableNames','Pressure');
%         fprintf('\n')
%         toc(pressure_timer);
%         fprintf('\n')
%     end

% Calculate the pressure using the time, voltage, and hydrophone senstivity curve
    function get_pressure
        pressure_timer=tic;
        fprintf('Converting voltages to pressure...\n')
        % Adjust the time to start at zero subtract the mean voltage (must be done after get_reference_indices)
        A.Time = A.Time-A.Time(:,1);
        A.Voltage = A.Voltage-mean(A.Voltage,2);
        
        [A.Time,A.Voltage]=zeropad_tv(A.Time,A.Voltage);
        
        p1=zeros(size(A.Time));
        fprintf('0%% -- ')
        calibration_data=readtable('HGL0200_calibration.csv');
%         calibration_data='HGL0200_calibration.csv';
        for n = 1:size(A.Time,1)
            p1(n,:) = reshape(calibrate_V2Pa(A.Time(n,:),A.Voltage(n,:),calibration_data),1,[]);
            
            % print progress
            myprogress = round(n/size(A.Time,1)*100);
            if any(n== round(size(A.Time,1)*(0:0.1:1)))
                fprintf('%03.0d%% -- ',myprogress)
            end
        end
        A=addvars(A, p1,'After','Time','NewVariableNames','Pressure');
        fprintf('\n')
        toc(pressure_timer);
        fprintf('\n')
    end


% Look at the two redundant measurements and use the one with the minimum rarefactive pressure as the reference
    function get_reference_indices
        ref_timer=tic;
        fprintf('Determining reference signals...\n')
        RefInds = zeros([size(A,1),1]);
        try
            for n=1:size(A,1)
                if A.Medium(n)=='rib' %#ok<BDSCA>
                    % identify possible reference signals
                    ids = A.Index((A.SampleNumber ==A.SampleNumber(n)) & (A.Medium=='saline') & (A.Power==A.Power(n)) & (A.frequency_MHz==A.frequency_MHz(n)) & (A.usProbe==A.usProbe(n)) & (A.date == A.date(n)) );
                    
                    % choose the reference signal with the lowest minimum voltage
                    voltsigs=A.Voltage(ids,:);
                    [~,mini]=min(min(voltsigs,[],2));
                    ids=ids(mini);
                    
                    RefInds(n)=ids;                                        
                    
                    % The negative of the reference number is assigned to the measurement with the higher PRPA, which will later be ignored.
                    ids = A.Index((A.SampleNumber ==A.SampleNumber(n)) & (A.Medium=='rib') & (A.Power==A.Power(n)) & (A.frequency_MHz==A.frequency_MHz(n)) & (A.usProbe==A.usProbe(n)) & (A.date == A.date(n)) );
                     % choose the reference signal with the lowest minimum voltage
                    voltsigs=A.Voltage(ids,:);
                    if min(A.Voltage(n,:)) ~= min(voltsigs(:))
                        RefInds(n) = -RefInds(n);
                    end
                end
                ids=[]; %#ok<NASGU>
                
                
                
                
            end
        catch myerror
            fprintf('Failed at n=%.0f',n)
            error(myerror)
        end
        A=addvars(A, RefInds,'After','tekNumber','NewVariableNames','ReferenceIndex');
        toc(ref_timer)
        fprintf('\n')
    end
end


