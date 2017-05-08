function l = GeneSubsetPower(m, GenesIn, SimMat, AllMu)
% l = GeneSubsetPower(m, Genes, SimMat)
% given a cluster model m, see how many hierarchical levels on average 
% the set Genes can distinguish
%
% Genes can be a numeric array, or a cell array of strings (gene names)
% SimMat is how similar the classes should be called

if isnumeric(GenesIn)
    Genes = GenesIn;
else
    
    Genes = zeros(length(GenesIn),1);

    for i=1:length(GenesIn)
        Genes(i) = find(strcmp(GenesIn{i}, m.GeneName));
    end
end
    

m2 = m;
% m2.Hard = 0;

m2.Active = Genes;
m2.mu = AllMu(Genes,:)';

m3 = m2.Estep();



%Sim = SimMat(sub2ind(size(SimMat), m.Class, m3.Class)); % nK by nK

% % figure(234); hist(Sim,0:6);
%l = mean(Sim);

%l = mean(sum(SimMat(m.Class,:).*m3.w',2));
w = LogLtoP(m3.L'); % prediction weights of new classes, nK by nC 
AvSim = w'*SimMat; % nK by nC 
l = mean(AvSim(sub2ind(size(AvSim),1:m.nC,m.Class(:)')));
