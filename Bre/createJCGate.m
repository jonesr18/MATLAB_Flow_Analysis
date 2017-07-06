function createJCGate(fileName)
% creates a file with gates for each of the Just Cell gates:
%       gate1: SSC-A vs. FSC-A
%       gate2: FSC-H vs. FSC-W
%       gate3: SSC-H vs. SSC-W


dirs = strfind(fileName,'/');
if dirs
    gateFileName = [fileName(1:dirs(end)) 'JCGate_' fileName(dirs(end)+1:end-4)];
else
    baseName=fileName(1:end-4);
    gateFileName=['JCGate_' baseName];
end

if exist([pwd '/' gateFileName '.mat'],'file')~=2

    [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(fileName);

    %P1
    ch_SSC_A=getChannel(fcshdr,'SSC-A');
    ch_FSC_A=getChannel(fcshdr,'FSC-A');

    datSSC_A=fcsdat(:,ch_SSC_A);
    datFSC_A=fcsdat(:,ch_FSC_A);

    figure
    xlabel('FSC-A')
    ylabel('SSC-A')
    hold on
    [P1inds gate1]=gatePolygon(datFSC_A,datSSC_A,'semilogy');
    

    %P2
    ch_FSC_W=getChannel(fcshdr,'FSC-W');
    ch_FSC_H=getChannel(fcshdr,'FSC-H');

    datFSC_W=fcsdat(P1inds,ch_FSC_W);
    datFSC_H=fcsdat(P1inds,ch_FSC_H);

    figure
    xlabel('FSC-W')
    ylabel('FSC-H')
    hold on
    [P2inds gate2]=gatePolygon(datFSC_W,datFSC_H,'linear');
    

    %P3
    ch_SSC_W=getChannel(fcshdr,'SSC-W');
    ch_SSC_H=getChannel(fcshdr,'SSC-H');

    datSSC_W=fcsdat(P1inds(P2inds),ch_SSC_W);
    datSSC_H=fcsdat(P1inds(P2inds),ch_SSC_H);

    figure
    xlabel('SSC-W')
    ylabel('SSC-H')
    hold on
    [P3inds gate3]=gatePolygon(datSSC_W,datSSC_H,'semilogy');
    

    JCInds=P1inds(P2inds(P3inds));
    
    save(gateFileName,'gate1','gate2','gate3')
else
%     warnmsg=[gateFileName '.mat already exists'];
%     warning(warnmsg)
end

end