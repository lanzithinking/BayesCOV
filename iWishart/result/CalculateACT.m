function [ ACTs ] = CalculateACT( Samples, Batchsize )


% Samples is a NumOfSamples x NumOfParameters matrix

[NumOfSamples, NumOfParameters] = size(Samples);

if(nargin<2) Batchsize = floor(NumOfSamples^(2/3)); end

NumOfBatches = floor(NumOfSamples/Batchsize);


% Calculate autocorrelation times
for i = 1:NumOfParameters
    ACTs(i) = act(Samples(NumOfSamples-Batchsize*NumOfBatches+1:end,i),Batchsize,NumOfBatches);
end


disp('ACT Values:')
disp(ACTs)


end


function [ACT] = act(Sample, Batchsize, NumOfBatches)

for b = 1:NumOfBatches
    Batchmeans(b) = mean(Sample(1+Batchsize*(b-1):Batchsize*b));
end

ACT = Batchsize*var(Batchmeans)/var(Sample);

end