%% SFDI Processing Code
% This script is designed to be run in sections. To run an individual
% section, simply click on the section to highlight it, then click the
% small "Run Section" button in the editor tab, or use (ctrl + enter).
% Read the code guide writup for a breif explaination of each code section.
% Each section should be run sequentially from top to bottom
%% Import data set for OP Calculations (from OpenSFDI system)
addpath('Chromophore_OPs','Data_Sets','Functions','LUTs','Phantom_OPs');
directory = uigetdir();
[dataSet,fileNames] = create_data_set(directory);
Names = string(fileNames);
for i = 1:length(Names)
    Names(i) = extractBefore(fileNames(i),'.');
    disp(strcat(num2str(i),': ',Names(i)));
end
%% Import data set for OP Calculations (from Modulim system)
addpath('Chromophore_OPs','Data_Sets','Functions','LUTs','Phantom_OPs');
fileIdx = 1; %This is the repetition ID (run this script for each repetition for some statistics)
directory = uigetdir();
[dataSet,fileNames] = create_data_set_MI(directory,fileIdx,[0,0.1],[0,900]);
Names = string(fileNames);
for i = 1:length(Names)
    disp(strcat(num2str(i),': ',Names(i)));
end
%% Load Calibration Phantom Data
phantomID = 1;
%%%%%%%%%%%%%
referencePhantom = dataSet{phantomID};
referencePhantom.Name = Names(phantomID);
referencePhantom = exposureCorrect(referencePhantom);
referencePhantom = demodulate(referencePhantom);
disp('Calibration Phantom Loaded');
%% Exposure Correct and Demodulate Raw Data
disp('Demodulating... ')
% Run this section to demodulate the raw data
progressbar('Demodulating... ');
for i = 1:length(dataSet)
    dataSet{i} = exposureCorrect(dataSet{i});
    dataSet{i} = demodulate(dataSet{i});
    progressbar(i/length(dataSet));
end
disp('Demodulation Finished');
%% Check Demodulation
measurementID = 2;
freqID = 2;
wavelengthID = 1;
%%%%%%%%%%%%%%%%
imagesc(dataSet{measurementID}.Demodulated(:,:,wavelengthID,freqID));
fx = dataSet{measurementID}.Freqs(freqID);
wv = dataSet{measurementID}.Wavelengths(wavelengthID);
title({'Demodulation Check',['f_x = ' num2str(fx) ', \lambda = ' num2str(wv)]});
%% Load Calibration File
opCalibFile = 'MI_Phantom_Optical_Properties';
ops = load([pwd '\Phantom_OPs\' opCalibFile]);
phantom_ops = ops.phantom_ops;
disp('Calibration File Loaded');
%% Calibrate Fix for Rd with LUT
LUTName = 'LUT_Gardner_homogeneous_0_0p1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\LUTs\' LUTName]);
disp('Calibrating... ')
Wavelengths = referencePhantom.Wavelengths;
Freqs = referencePhantom.Freqs;
ACph = referencePhantom.Demodulated;
progressbar('Calibrating... ');
for i = 1:length(dataSet)
    dataSet{i}.Rd = calibrate_new(Wavelengths,Freqs,dataSet{i}.Demodulated,ACph,LUT,phantom_ops.data);
    progressbar(i/length(dataSet));
end
disp('Rd Calculated');
%% Calibrate with Diffusion Model
opCalibFile = 'MI_Phantom_Optical_Properties';
ops = load([pwd '\Phantom_OPs\' opCalibFile]);
phantom_ops = ops.phantom_ops; phantom_ops = phantom_ops.data;
disp('Calibration File Loaded');
for i = 1:length(dataSet)
    for j = 1:length(dataSet{i}.Wavelengths)
        muaPhantom = interp1(phantom_ops(:,1), phantom_ops(:,2),dataSet{i}.Wavelengths(j),'pchip', 'extrap');
        muspPhantom = interp1(phantom_ops(:,1), phantom_ops(:,3),dataSet{i}.Wavelengths(j),'pchip', 'extrap');
        for k = 1:length(dataSet{i}.Freqs)
            RdModel = diffApprox(muaPhantom,muspPhantom,1.4,dataSet{i}.Freqs(k));
            dataSet{i}.Rd(:,:,j,k) = (dataSet{i}.Demodulated(:,:,j,k)./dataSet{phantomID}.Demodulated(:,:,j,k))*RdModel;
        end
    end
    progressbar(i/length(dataSet));
end
%% Plot Rd-Freq Relation
measurementID = 1;
wavelengthID = 3;
yroi = 320:340;
xroi = 250:270;
Rd = squeeze(dataSet{measurementID}.Rd(yroi,xroi,wavelengthID,:));
RdMean = squeeze(mean(mean(Rd)));
plot(dataSet{measurementID}.Freqs,RdMean,'-o','LineWidth',2);
hold on
%% Height Correction
heightFactor = -2.0975; %(1/hP(1))
RdCorrectionFactors = [0.0095,0.0025;0.0051,0.0009;0.0088,0.0017];
%RdCorrectionFactors = squeeze(rdP(1,:,:)); %rdP(slope,wv,fx)
pFreqID = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(dataSet)
    dataSet{i}.Height = heightFactor*getProfile(dataSet{phantomID}.Raw.Profile(:,:,pFreqID,:),dataSet{i}.Raw.Profile(:,:,pFreqID,:));
    for j = 1:length(dataSet{phantomID}.Freqs)
        for k = 1:length(dataSet{phantomID}.Wavelengths)
            dataSet{i}.RdHC(:,:,k,j) = dataSet{i}.Rd(:,:,k,j) - dataSet{i}.Height*RdCorrectionFactors(k,j);
        end
    end
end
disp('Height Correction Finished');
%% Plot Height Map
measurementID = 1;
yroi = 300:800;
xroi = 300:900;
%%%%%%%%%%%%%%%
figure
mesh(dataSet{phantomID}.Height(yroi,xroi));
hold on
mesh(dataSet{measurementID}.Height(yroi,xroi));
title('Height Map of Sample Relative to Phantom');
zlabel('Height (mm)');
%% Plot Corrected Rd
measurementID = 1;
freqID = 1;
wavelengthID = 1;
yroi = 300:800;
xroi = 300:900;
%%%%%%%%%%%%%%%
figure
imagesc(100*abs(dataSet{measurementID}.RdHC(yroi,xroi,wavelengthID,freqID)-dataSet{measurementID}.Rd(yroi,xroi,wavelengthID,freqID))./dataSet{measurementID}.Rd(yroi,xroi,wavelengthID,freqID));
title('Rd Corrected vs. Non-Corrected');
cb = colorbar;
title(cb,'% Difference');
%% OP Lookup
LUTName = 'LUT_Gardner_homogeneous_0_0p1';
load(['LUTs\',LUTName]);
%If height corrected use 1, if not use 0
useCorrectedRd = 0;
yroi = 1:696;
xroi = 1:520;
%%%%%%%%%%%%%%%%%%%
if useCorrectedRd == 1
    disp('Calculating OPs... ')
operations = length(dataSet)*length(Wavelengths);
count = 0;
progressbar('Calculating Optical Properties of Data Set...');
for i = 1:length(dataSet)
    for j = 1:length(Wavelengths)
        dataSet{i}.OP(:,:,j,1) = griddata(LUT.M1,LUT.M2,LUT.Mua,dataSet{i}.RdHC(yroi,xroi,j,1),dataSet{i}.RdHC(yroi,xroi,j,2));
        dataSet{i}.OP(:,:,j,2) = griddata(LUT.M1,LUT.M2,LUT.Musp,dataSet{i}.RdHC(yroi,xroi,j,1),dataSet{i}.RdHC(yroi,xroi,j,2));
        count = count + 1;
        progressbar(count/operations);
    end
end
else
disp('Calculating OPs... ')
operations = length(dataSet)*length(Wavelengths);
count = 0;
progressbar('Calculating Optical Properties of Data Set...');
for i = 1:length(dataSet)
    for j = 1:length(Wavelengths)
        dataSet{i}.OP(:,:,j,1) = griddata(LUT.M1,LUT.M2,LUT.Mua,dataSet{i}.Rd(yroi,xroi,j,1),dataSet{i}.Rd(yroi,xroi,j,2));
        dataSet{i}.OP(:,:,j,2) = griddata(LUT.M1,LUT.M2,LUT.Musp,dataSet{i}.Rd(yroi,xroi,j,1),dataSet{i}.Rd(yroi,xroi,j,2));
        count = count + 1;
        progressbar(count/operations);
    end
end
end
disp('Optical Properties Calculated');
%% Calculate A and B maps w/ error
k = 1:length(dataSet); k(phantomID) = [];
for k = k
    [dataSet{k}.a,dataSet{k}.b,dataSet{k}.abError] = getAandBMaps(dataSet{k});
    progressbar(k/(length(dataSet)-1));
end
%% Hemoglobin Calculations
M = chromophore_op_lookup(dataSet{phantomID}.Wavelengths);
HbO = M(:,1)';
Hb = M(:,2)';
[y,x] = size(dataSet{phantomID}.OP(:,:,1,1));
N = length(dataSet)-1;
k = 1:length(dataSet); k(phantomID) = [];
c = 0;
for k = k
c = c+1;
disp(['working on ',num2str(c),' of ',num2str(N)]);
[y,x,n] = size(dataSet{phantomID}.OP(:,:,1,1));
dataSet{k}.OxyMap = zeros(y,x);
dataSet{k}.deOxyMap = zeros(y,x);
dataSet{k}.OxyDeoxyError = zeros(y,x);
for i = 1:y
    for j = 1:x
        absorption = squeeze(dataSet{k}.OP(i,j,wvIdx,1))';
        [alpha,beta,error] = hemoglobinFit(HbO,Hb,absorption);
        dataSet{k}.OxyMap(i,j) = alpha;
        dataSet{k}.deOxyMap(i,j) = beta;
        dataSet{k}.OxyDeoxyError(i,j) = error;
    end
    progressbar(i/y);
end
end
disp('Hemoglobin Concentrations Calculated');
%% Save Processed Data Set
dataSetName = 'Rosemary_Hand_ROI_LUT_Method';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fields = {'Raw','Phases'};
for i = 1:length(dataSet)
    dataSet{i} = rmfield(dataSet{i},fields);
end
save(dataSetName,'dataSet','-v7.3');
disp([dataSetName 'saved to' pwd]);