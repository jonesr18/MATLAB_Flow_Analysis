function datYnew_final = comp_mixed(sample,chbleed,ch2corr,scBleed,sJ,plotting)
% compensates using amixed model based on logicle transform

if ~exist('plotting')
    plotting = 0;
end

% calculate autoflourescence values and transfection cutoffs using
% nontransfected sample
    [dat0, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(sJ);
    chY=ch2corr;%getChannel(fcshdr,'FIT');
    chR=chbleed;%getChannel(fcshdr,'Red');
    chB=getChannel(fcshdr,'Blue');
    datY = dat0(:,chY);
    datR = dat0(:,chR);
    datB = dat0(:,chB);
    autoF_Y = mean(datY);
    autoF_R = mean(datR);
    autoF_B = mean(datB);
    txY = prctile(datY-autoF_Y,99.9);
    txR = prctile(datR-autoF_R,99.9);
    txB = prctile(datB-autoF_B,99.9);

%SC
%--> Red bleeds into yellow

% solve for mf using single color sample
% calculate mf where positively transfected cells in bleedthrough color
% have population centered at zero in channel correcting in

    %subtract autoflourescence
    dat = applyJCGate(scBleed,sJ);
    datY = dat(:,chY)-autoF_Y;
    datR = dat(:,chR)-autoF_R;
    datB = dat(:,chB)-autoF_B;
    
    if plotting
    figure
    plot(lin2logicle(datY),lin2logicle(datR),'.');
    biexaxis(gca)
    end
    
    
    mf = fsolve(@minimizer,0.009);
    mf = 0.367;
    
    function out = minimizer(mf)
        datYn = datY - mf.*datR;
        out = mean(datYn(datR>txR*10));
    end

    datYnew = datY - mf.*datR;
                                                            if plotting
                                                                figure
                                                                plot(lin2logicle(datYnew),lin2logicle(datR),'.');
                                                                biexaxis(gca)
                                                            end
        

% for given sample:

    dat = applyJCGate(sample,sJ);
    datY = dat(:,chY)-autoF_Y;
    datR = dat(:,chR)-autoF_R;
    datB = dat(:,chB)-autoF_B;

                                                            if plotting
                                                                %original data
                                                                figure
                                                                subplot(4,3,1)
                                                                plot(lin2logicle(datR),lin2logicle(datY),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                                xlabel('Red')
                                                                ylabel('Yellow')
                                                                subplot(4,3,2)
                                                                plot(lin2logicle(datB),lin2logicle(datY),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                                xlabel('Blue')
                                                                ylabel('Yellow')
                                                                subplot(4,3,3)
                                                                plot(lin2logicle(datR),lin2logicle(datB),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                                xlabel('Red')
                                                                ylabel('Blue')
                                                            end
        
        datYnew = datY - datR.*mf;
        
                                                            if plotting
                                                    %         figure
                                                                subplot(4,3,4)
                                                                plot(lin2logicle(datR),lin2logicle(datYnew),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                                subplot(4,3,5)
                                                                plot(lin2logicle(datB),lin2logicle(datYnew),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                                subplot(4,3,6)
                                                                plot(lin2logicle(datR),lin2logicle(datB),'.')
                                                                axis([0 4.5 0 4.5])
                                                                biexaxis(gca)
                                                            end
        
        
        datYnew = logicle2lin(lin2logicle(datY) - lin2logicle(datR.*mf) + lin2logicle(autoF_Y));

        if plotting
%         figure
        subplot(4,3,7)
        plot(lin2logicle(datR),lin2logicle(datYnew),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
        subplot(4,3,8)
        plot(lin2logicle(datB),lin2logicle(datYnew),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
        subplot(4,3,9)
        plot(lin2logicle(datR),lin2logicle(datB),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
        end
        
        
        datY_LN = datY - datR.*mf;
        datY_LG = logicle2lin(lin2logicle(datY) - lin2logicle(datR.*mf) + lin2logicle(autoF_Y));
        T=262144;   % The top data value
        M=4.5;      % Breadth of the display in decades
        r=-20;%-150;    % Negative range reference value
        W=(M-log10(T/abs(r)))/2
%         syms p
%         P=solve(W==2*p*log10(p)/(p+1),p)
        P = 2.6973738273810843247742159802084;
        X = (lin2logicle(datY) - lin2logicle(datR.*mf)).*(4.5 - lin2logicle(autoF_Y))./(4.5 - lin2logicle(datR.*mf)) + lin2logicle(autoF_Y);
%         X = XD(XD>=W);
%         datY_LGu = datY_LG(XD>=W);
%         datY_LNu = datY_LN(XD>=W);
%         datRu = datR(XD>=W);
%         datBu = datB(XD>=W);
        datYnew = datY_LN.*10.^(X-W)./(10.^(X-W) + P^2.*10.^(-(X-W)./P)) + datY_LG.*P^2.*10.^(-(X-W)./P)./(10.^(X-W) + P^2.*10.^(-(X-W)./P));
%         datRu = datR(XD>=W);
%         datBu = datB(XD>=W);
%         datYnewu = datY_LN.*10.^(X-W)./(10.^(X-W) + P^2.*10.^(-(X-W)./P)) + datY_LG.*P^2.*10.^(-(X-W)./P)./(10.^(X-W) + P^2.*10.^(-(X-W)./P));
        
if plotting
        subplot(4,3,10)
        plot(lin2logicle(datR),lin2logicle(datYnew),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
        subplot(4,3,11)
        plot(lin2logicle(datB),lin2logicle(datYnew),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
        subplot(4,3,12)
        plot(lin2logicle(datR),lin2logicle(datB),'.')
        axis([0 4.5 0 4.5])
        biexaxis(gca)
end

        datYnew_final = datYnew + autoF_Y;

end