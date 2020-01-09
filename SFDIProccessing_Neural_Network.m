%% SFDI Processing
%% Import data set for OP Calculations (from OpenSFDI system)
addpath('Chromophore_OPs','Data_Sets','Functions','LUTs','Phantom_OPs','NNModels');
directory = uigetdir();
[dataSet,fileNames] = create_data_set(directory);
Names = string(fileNames);
for i = 1:length(Names)
    Names(i) = extractBefore(fileNames(i),'.');
    disp(strcat(num2str(i),': ',Names(i)));
end
%% Import data set for OP Calculations (from Modulim system)
addpath('Chromophore_OPs','Data_Sets','Functions','LUTs','Phantom_OPs','NNModels');
directory = uigetdir();
fileIdx = 1; %This is the repetition ID (run this script for each repetition for some statistics)
[dataSet,fileNames] = create_data_set_MI(directory,fileIdx,[0,0.05,0.1,0.15,0.2],[1,1000]);
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
measurementID = 3;
freqID = 2;
wavelengthID = 8;
%%%%%%%%%%%%%%%%
imagesc(dataSet{measurementID}.Demodulated(:,:,wavelengthID,freqID));
fx = dataSet{measurementID}.Freqs(freqID);
wv = dataSet{measurementID}.Wavelengths(wavelengthID);
title({'Demodulation Check',['f_x = ' num2str(fx) ', \lambda = ' num2str(wv)]});
%% Calibrate with Forward Model
opCalibFile = 'MI_Phantom_Optical_Properties';
ops = load([pwd '\Phantom_OPs\' opCalibFile]);
phantom_ops = ops.phantom_ops; phantom_ops = phantom_ops.data;
disp('Calibration File Loaded');
for i = 1:length(dataSet)
    for j = 1:length(dataSet{i}.Wavelengths)
        muaPhantom = interp1(phantom_ops(:,1), phantom_ops(:,2),dataSet{i}.Wavelengths(j),'pchip', 'extrap');
        muspPhantom = interp1(phantom_ops(:,1), phantom_ops(:,3),dataSet{i}.Wavelengths(j),'pchip', 'extrap');
        for k = 1:length(dataSet{i}.Freqs)
            RdModel = getRdHomogeneous([muaPhantom,muspPhantom],dataSet{i}.Freqs(k));
            dataSet{i}.Rd(:,:,j,k) = (dataSet{i}.Demodulated(:,:,j,k)./dataSet{phantomID}.Demodulated(:,:,j,k))*RdModel;
        end
    end
    progressbar(i/length(dataSet));
end
%% OP Calculation
yroi = 240:420;
xroi = 180:340;
dim = [length(yroi),length(xroi)];
for i = 1:length(dataSet)
    for j = 1:length(dataSet{i}.Wavelengths)
    Rd = zeros(numel(dataSet{i}.Rd(yroi,xroi,j,1)),size(dataSet{i}.Rd(yroi,xroi,j,:),4));
        for k = 1:length(dataSet{i}.Freqs)
            Rd(:,k) = reshape(dataSet{i}.Rd(yroi,xroi,j,k),1,[]);
        end
        OPs = HomogeneousInverse(Rd);
        dataSet{i}.OP(:,:,j,1) = reshape(OPs(:,1),dim);
        dataSet{i}.OP(:,:,j,2) = reshape(OPs(:,2),dim);
    end
    progressbar(i/length(dataSet));
end
%% Calculate A and B maps w/ error
yroi = 310:380;
xroi = 220:290;
for i = 1:length(dataSet)
    [dataSet{i}.a,dataSet{i}.b,dataSet{i}.abError] = getAandBMaps(dataSet{i});
    progressbar(i/length(dataSet));
end
%% Hemoglobin Calculations
M = chromophore_op_lookup(dataSet{phantomID});
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
        absorption = squeeze(dataSet{k}.OP(i,j,:,1))';
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
dataSetName = 'Rosemary_Hand_ROI_NN_Method';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fields = {'Raw','Phases'};
for i = 1:length(dataSet)
    dataSet{i} = rmfield(dataSet{i},fields);
end
save(dataSetName,'dataSet','-v7.3');
disp([dataSetName ' saved to: ' pwd]);