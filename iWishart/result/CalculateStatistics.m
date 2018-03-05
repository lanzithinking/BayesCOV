function [ESS,ACT,Times,Acpts] = CalculateStatistics(filabl,resloc)

% Get stats
Files = dir([resloc,'*' filabl '*.mat']);


for i = 1:length(Files)

    Data = open(Files(i).name);
    
    ESS.raw{i} = CalculateESS(Data.samp, size(Data.samp,1)-1);
    
    ESS.min(i) = min(ESS.raw{i});
    ESS.max(i) = max(ESS.raw{i});
    ESS.med(i) = median(ESS.raw{i});
    ESS.mean(i) = mean(ESS.raw{i});
    
    ACT.raw{i} = CalculateACT(Data.samp);
    ACT.min(i) = min(ACT.raw{i});
    ACT.max(i) = max(ACT.raw{i});
    ACT.med(i) = median(ACT.raw{i});
    ACT.mean(i) = mean(ACT.raw{i});
    
    Times(i) = Data.time;
    Acpts(i) = Data.acpt(1);
    
end

disp(['Time:   ' num2str(mean(Times)) ' +/- ' num2str(std(Times)/sqrt(length(Times)))])
disp(' ')

disp(['ESS for ' filabl ' dataset.'])
disp(' ')

disp(['Min:    ' num2str(mean(ESS.min)) ' +/- ' num2str(std(ESS.min)/sqrt(length(ESS.min)))])
disp(['Median: ' num2str(mean(ESS.med)) ' +/- ' num2str(std(ESS.med)/sqrt(length(ESS.med)))])
disp(['Mean:   ' num2str(mean(ESS.mean)) ' +/- ' num2str(std(ESS.mean)/sqrt(length(ESS.mean)))])
disp(['Max:    ' num2str(mean(ESS.max)) ' +/- ' num2str(std(ESS.max)/sqrt(length(ESS.max)))])

disp('')
disp(['Min ESS per seconds: ' num2str(mean(ESS.min)/mean(Times))])
disp(' ')


disp(['ACT for ' filabl ' dataset.'])
disp(' ')

disp(['Min:    ' num2str(mean(ACT.min)) ' +/- ' num2str(std(ACT.min)/sqrt(length(ACT.min)))])
disp(['Median: ' num2str(mean(ACT.med)) ' +/- ' num2str(std(ACT.med)/sqrt(length(ACT.med)))])
disp(['Mean:   ' num2str(mean(ACT.mean)) ' +/- ' num2str(std(ACT.mean)/sqrt(length(ACT.mean)))])
disp(['Max:    ' num2str(mean(ACT.max)) ' +/- ' num2str(std(ACT.max)/sqrt(length(ACT.max)))])

disp('')
disp(['Max ACT multiply time per iteration: ' num2str(mean(ACT.max)*mean(Times)/size(Data.samp,1))])

end
