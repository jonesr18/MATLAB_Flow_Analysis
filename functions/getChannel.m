function ch=getChannel(fcshdr,strName)
    % GETCHANNEL(FCSHDR,STRNAME) gets the channel number from fcs data
    %   CH = GETCHANNEL(FCSHDR, STRNAME) returns the column number
    %   corresponding to the color channel indicated by strName, where
    %    FCSHDR - Header data obtained from FCA_READFCS
    %    STRNAME - String used to uniquely identify the channel name
    %    CH - The index of the given channel
    %
    %   Example:
    %       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
    %       GreenChannel = getChannel(fcshdr,'FIT');
    %       GreenData = fcsdat(:,GreenChannel);
    %
    %   Written by
    %   Breanna Stillo
    %   bstillo@mit.edu
    %   Last Updated: 2014-10-14

    ch = find(strcmpi({fcshdr.par.name}, strName));
    if (isempty(ch))
        error(['chanel for "' strName '" not found. Try a different string identifier'])
    end
end