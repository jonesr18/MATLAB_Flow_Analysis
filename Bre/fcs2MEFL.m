function out = fcs2MEFL(fcsdat,channelNum,beadCalibrationFile)
%out = fcs2MEFL(fcsdat,channelNum,beadCalibrationFile)
%converts fcs data to MEFL
%   fcsdat = the fcsdat to be converted (may be a single column vector or a
%           matrix where each column corresponds to a different channel
%   channelNum = the channel numbers corresponding to the fcsdat columns.
%           Must be specified if the number of channels in fscdat is less than that
%           in the calibration file
%   beadCalibrationFile = the bead calibration file obtained from running
%           'calibrateBeads.m'.  If this is not specified, the user can
%           select the file from a menu
%   out = output data the same size as fcsdat, corresponding to MEFL
%           converted values
%
%   Example:
%       calibrateBeads
%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
%       GreenChannel = getChannel(fcshdr,'FIT');
%       RedChannel = getChannel(fcshdr,'Red');
%       GreenData = fcsdat(:,GreenChannel);
%       RedData = fcsdat(:,RedChannel);
%       MEFLdata = fcs2MEFL([GreenData RedData],[GreenChannel RedChannel])
%       GreenMEFL = MEFLdata(:,1);
%       RedMEFL = MEFLdata(:,2);
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-22;

switch nargin
    % only fcsdat is supplied
    case 1
        beadCalibrationFile = uigetfile('*.mat','Select the Beads Calibration File for this experiement');
        varName=beadCalibrationFile(1:end-4);
        fitVar=load([varName '/' beadCalibrationFile]);
        ChannelFits=fitVar.(varName);
        NumTotalChannels=length(ChannelFits(:,1));
        NumDataChannels=length(fcsdat(1,:));
        if NumDataChannels==1
            error('Channel number must be supplied')
        elseif NumDataChannels~=NumTotalChannels
            error('The number of channels in the dataset does not match thr number of channels in the calibration file.  Either supply the channel number of the given data set or supply a data set for all channels')
        else
            channelNum=1:length(fcsdat(1,:));
        end
    % fcsdat and channelNum are supplied which is either the list of channel numbers or the bead
    % calibration file
    case 2
        % check whether string (calibration file) or channel numbers
        if strcmp(class(channelNum),'char') % it is the calibration file
            beadCalibrationFile=channelNum;
            varName=beadCalibrationFile(1:end-4);
            fitVar=load([varName '/' beadCalibrationFile]);
            ChannelFits=fitVar.(varName);
            NumTotalChannels=length(ChannelFits(:,1));
            NumDataChannels=length(fcsdat(1,:));
            if NumDataChannels==1
                error('Channel number must be supplied')
            elseif NumDataChannels~=NumTotalChannels
                error('The number of channels in the dataset does not match thr number of channels in the calibration file.  Either supply the channel number of the given data set or supply a data set for all channels')
            else
                channelNum=1:length(fcsdat(1,:));
            end
        else % it was the channel numbers so ask user for calibration file
            beadCalibrationFile = uigetfile('*.mat','Select the Beads Calibration File for this experiement');
            varName=beadCalibrationFile(1:end-4);
            fitVar=load([varName '/' beadCalibrationFile]);
            ChannelFits=fitVar.(varName);
            NumTotalChannels=length(ChannelFits(:,1));
            NumDataChannels=length(fcsdat(1,:));
        end
    % all arguments are supplied
    case 3
        if strcmp(class(beadCalibrationFile),'char') % it is the calibration file
            varName=beadCalibrationFile(1:end-4);
            fitVar=load(beadCalibrationFile);
            ChannelFits=fitVar.(varName);
            NumTotalChannels=length(ChannelFits(:,1));
            NumDataChannels=length(fcsdat(1,:));
        else % it is the matrix
            ChannelFits = beadCalibrationFile;
            NumTotalChannels=length(ChannelFits(:,1));
            NumDataChannels=length(fcsdat(1,:));
        end
end

r=2^18; %resolution
n=4.5; %log decades


%test for negative values
ind=logical(sum(fcsdat<=0,2));
NumNegs=sum(ind);

if NumNegs>0
    msg=['Attempted to convert negative values to MEFL. ' num2str(NumNegs) ' of ' num2str(length(fcsdat(:,1))) ' cells were removed.'];
    warning(msg)
    fcsdat=fcsdat(~ind,:);
end

newdata=zeros(size(fcsdat));

ind=1;
for ch=channelNum
    ch
    intensity=fcsdat(:,ind);
    alpha=ChannelFits(ch,1);
    beta=ChannelFits(ch,2);
%     newdata(:,ind)=10^beta.*intensity.^(alpha*r/n); %convert to MEFL
    newdata(:,ind) = 10.^polyval(ChannelFits(ch,:),log10(intensity)); %convert to MEFL
    ind=ind+1;
end
out=newdata;

end