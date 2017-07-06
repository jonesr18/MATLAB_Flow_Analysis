function fcsdatGated = applyJCGate(sampleFile, gateFile)
% applies a certain JC gate to a specified sample
    suffix = gateFile(end-2:end);
    if strcmp(suffix, 'fcs')
        createJCGate(gateFile);
        
        dirs = strfind(gateFile,'/');
        if dirs
            gateFileName = [gateFile(1:dirs(end)) 'JCGate_' gateFile(dirs(end)+1:end-4) '.mat'];
        else
            baseName=gateFile(1:end-4);
            gateFileName=['JCGate_' baseName '.mat'];
        end
        
%         baseName=gateFile(1:end-4);
%         gateFileName=['JCGate_' baseName '.mat'];
    else
        gateFileName=gateFile;
    end
        
    GF=load(gateFileName);

    [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(sampleFile);
    ch_SSC_A=getChannel(fcshdr,'SSC-A');
    ch_FSC_A=getChannel(fcshdr,'FSC-A');
    ch_FSC_W=getChannel(fcshdr,'FSC-W');
    ch_FSC_H=getChannel(fcshdr,'FSC-H');
    ch_SSC_W=getChannel(fcshdr,'SSC-W');
    ch_SSC_H=getChannel(fcshdr,'SSC-H');

    [P1inds g1]=gatePolygon(fcsdat(:,ch_FSC_A),fcsdat(:,ch_SSC_A),'semilogy',GF.gate1);
    [P2inds g2]=gatePolygon(fcsdat(P1inds,ch_FSC_W),fcsdat(P1inds,ch_FSC_H),'linear',GF.gate2);
    [P3inds g3]=gatePolygon(fcsdat(P1inds(P2inds),ch_SSC_W),fcsdat(P1inds(P2inds),ch_SSC_H),'semilogy',GF.gate3);

    JCInds=P1inds(P2inds(P3inds));
    
    fcsdatGated=fcsdat(JCInds,:);

end