function orthogonalityMatrix(A,labels)
%orthogonalityMatrix(A,labels)
%generates an orthogonality plot for the matrix A
%   A = matrix of data vlaues (must be square)
%   labels = labels for the matrix
%
%Ex:
%   A =[ 4.5000    0.5000    0.2500    0.5500
%     0.0500    3.5000    1.0000    0.2500
%     0.5000    0.7500    4.2500    0.1500
%     0.6500    0.4500    1.2500    4.7500];
%   orthoMatrix(A,{'A1','A2','A3','A4'})
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14;


map=[linspace(0,1)' zeros(100,1) linspace(1,0)'];
colormap(map)
imagesc(A)
colorbar('southoutside')

set(gca,'Xtick',1:length(A),...
    'Ytick',1:length(A),...
    'XAxisLocation','top')

    if exist('labels','var')
        set(gca,'XTickLabel',labels,'YTickLabel',labels)
    end
end