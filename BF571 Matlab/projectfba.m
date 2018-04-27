% renal.c(3913) = 1;
% Liver.c(4386) = 1;
% Lung.c(3809) = 1;
% Urothelial.c(4205) = 1;
Breast.c(4299) = 1;

model = Breast;
% optimize the biomass

% get a default analysis w/ Segre's FBA
R = FBA2(model.S,model.lb,model.ub,model.c,model.mets,model.rxns) ;
%declare empty cell that is preallocated for speed.  Size is arbirtrary
delGenes = cell(50,2);
% declare counter
c = 1;
for i = 1:length(model.genes)
    % use removeGenes function from Raven Toolbox
    [tempM,notdel] = removeGenes(model,i);
    %Use Segre's function 
    [T T_o T_s] = FBA2(tempM.S,tempM.lb,tempM.ub,tempM.c,tempM.mets,tempM.rxns) ;
    if T_o < .1
        % add any gene knockout to the matrix that satisfies the threshold
        delGenes(c,1) = {i};
        delGenes(c,2) = model.genes(i);
        % increase counter for gene cell structure
        c = c+1;
    end
    %percentage tracker 
    per = i*100/length(model.genes);
    fprintf('%.2f percent complete \n',per);
end
fprintf('count is: %d', c -1)


