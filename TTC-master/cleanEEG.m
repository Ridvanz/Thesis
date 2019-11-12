function [input,output] = cleanEEG(EEGdata,partici)

ds = 8;

PAR = size(EEGdata,1);          % #participants 
REA = size(EEGdata{2}.angle,1); % #realizations


for j = 1:REA
    out = squeeze(EEGdata{partici}.comp(j,:,:))';
    in  = squeeze(EEGdata{partici}.angle(j,:,:))';

outlierz = isoutlier(mean(out)) | isoutlier(var(out));

out=smoothdata(out(:,~outlierz),'gaussian',4*ds);
in=smoothdata(in(:,~outlierz),'gaussian',4*ds);

outbounds = (max(abs(out))<4);

EEGclean{partici}.out{j}=out(1:ds:end,outbounds)/8+0.5;
EEGclean{partici}.in{j}=in(1:ds:end,outbounds)/0.12+0.5;
    end

input = [];
output = [];
    for j = 1:REA
input = [EEGclean{partici}.in{j} input];
output = [EEGclean{partici}.out{j} output];
    end


end

