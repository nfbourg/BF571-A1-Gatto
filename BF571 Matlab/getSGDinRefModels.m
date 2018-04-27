function [toxicDataset,modelOpt,geneEssentiality,levyReport] = getSGDinRefModels(refModelName,medium,wFluxes,cellLine)   
    if nargin < 4 || isempty(cellLine)
        cellLine = '';
    end    
    perturbation       = 'sgd';
    switch medium
        case 'FBS'
            mediumMets = getOldFBSMets;
        case 'HAM'            
            mediumMets = getHAMMets;
    end
    parameters.suppressOutput      = false;
    parameters.useControl          = false;
    parameters.fixUnfeasibleFluxes = false;  
    parameters.convertBoundsToInfinite = false;
    parameters.genePerturbationMethod  = 'fba';
    parameters.saveWTmodel             = strcat('Models/Constrained/',refModelName,wFluxes,'_',cellLine,'_',medium,'.mat');
    parameters.saveSGD                 = strcat('SGDpredictions/SGDresultFiles/','SGDin',refModelName,wFluxes,'_',cellLine,'_',medium,'.txt');
    toxicThreshold = 0.5;
    switch wFluxes
        case 'woFluxes'
        objBounds = [0,1]; % -Inf yields far more toxic KOs
        typical   = false;
        switch refModelName
            case 'HMR3674'
                model = [];
                load HMR3674_20130326.mat
                model = addDiffusionCORErxns(model);
                model = convertBoundsToFinite(model,1000);
                obj  = {'HMAR_BiomassOut'};
                probMets = {'m01451s'
                            'm01830c'
                            'm02348c'
                            'm02402c'
                            'm02440c'
                            'm02583c'
                            'm02586c'
                            'm02588c'
                            'm02646s'
                            'm02783c'
                            'm02949c'
                            'm02965c'
                            'm02969c'
                            'm02998c'};
                model   = addExchangeRxns(model,'in',getIndexes(model,probMets,'mets'));
                parameters.singleGenesNotToPerturb = {};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality]  = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd')
            case 'test'
                load Models/matlab.mat
                model = addDiffusionCORErxns(renal);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CancerBiomass_Ex'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality]  = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
            case 'iRC1410'
                load Models/iRenalCancer1410.mat
                model = addDiffusionCORErxns(RCCmodel);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CancerBiomass_Ex'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality]  = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
            case 'test'
                model = []
                load Models/matlab.mat
                model = addDiffusionCORErxns(model);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CancerBiomass_Ex'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality]  = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
            case 'iPC1675'
                load Models/iProstateCancer1675.mat
                model = addDiffusionCORErxns(PROSmodel);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CancerBiomass_Ex'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
            case '786Ov1'
                load Models/786O_v1
                model = addDiffusionCORErxns(VHLmodel);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CANCERBOUT'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd')
            case '786Ov2'
                load Models/786O_v2
                model = addDiffusionCORErxns(VHLmodel);
                model = convertBoundsToFinite(model,1000);
                obj   = {'CANCERBOUT'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd')
            case 'PC3'
                load Models/cPC3
                model = addDiffusionCORErxns(PC3);
                model = convertBoundsToFinite(model,1000);
                model = addCompartmentalizedMetNames(model);
                mediumMets = unique([mediumMets(ismember(mediumMets,model.metNamesC));'linolenate[s]']);
                obj   = {'CANCERBOUT'};
                [modelOpt,geneEssentiality,levyReport] = assessModelViability(model,obj,perturbation,mediumMets,repmat([-1000 1000],numel(mediumMets),1),objBounds,typical,parameters);
                [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',medium));
                drawPerturbationSummary(geneEssentiality,'sgd')
        end
        case 'wFluxes'
            load COREforLeviathan.mat
            typical = false;
            switch refModelName
                case 'HMR3674'
                    load COREforLeviathan
                    load Models/HMR3674_20130326.mat
                    probMets = {'m01451s'
                                'm01830c'
                                'm02348c'
                                'm02402c'
                                'm02440c'
                                'm02583c'
                                'm02586c'
                                'm02588c'
                                'm02646s'
                                'm02783c'
                                'm02949c'
                                'm02965c'
                                'm02969c'
                                'm02998c'};
                    refModel     = addDiffusionCORErxns(model);
                    refModel     = addExchangeRxns(refModel,'in',getIndexes(refModel,probMets,'mets'));
                    refModel     = convertBoundsToFinite(refModel,1000);
                    refModel     = fixMisannotatedGene(refModel);
                    RCCcellLines = unique(CLs(strncmp('Renal',Histology,5)));
                    obj          = {'HMAR_BiomassOut'};
                    [~,allMetFluxes,metList]  = retrieveCommonExchangeFluxes(refModel,RCCcellLines,2,obj,true,true,medium);
                    metFluxes                 = allMetFluxes(:,:,ismember(RCCcellLines,cellLine)); 
                    objBounds                 = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured      = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);
                    [toxicDataset,geneEssentiality]   = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes  = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
                    filename = strcat('SGDpredictions/',cellLine,'_in',refModelName,'_',medium);
                    save(filename,'modelOpt','geneEssentiality','levyReport')
                case 'iRC1410'
                    load Models/iRenalCancer1410.mat
                    load COREforLeviathan
                    refModel     = addDiffusionCORErxns(RCCmodel);
                    refModel     = convertBoundsToFinite(refModel,1000);
                    refModel     = fixMisannotatedGene(refModel);
                    RCCcellLines = unique(CLs(strncmp('Renal',Histology,5)));
                    obj          = {'CancerBiomass_Ex'};
                    [~,allMetFluxes,metList]  = retrieveCommonExchangeFluxes(refModel,RCCcellLines,2,obj,true,true,medium);
                    metFluxes                 = allMetFluxes(:,:,ismember(RCCcellLines,cellLine)); 
                    objBounds                 = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured      = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);
                    [toxicDataset,geneEssentiality]   = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes  = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
                    filename = strcat('SGDpredictions/',cellLine,'_inRCCpanel_',medium);
                    save(filename,'modelOpt','geneEssentiality','levyReport')
                case 'iPC1675'
                    load Models/iProstateCancer1675.mat
                    load COREforLeviathan
                    refModel     = addDiffusionCORErxns(PROSmodel);
                    refModel     = convertBoundsToFinite(refModel,1000);
                    refModel     = fixMisannotatedGene(refModel);
                    PCcellLines  = unique(CLs(strncmp('Prostate',Histology,5)));
                    obj          = {'CancerBiomass_Ex'};
                    [~,allMetFluxes,metList]  = retrieveCommonExchangeFluxes(refModel,PCcellLines,2,obj,true,true,medium);
                    metFluxes                 = allMetFluxes(:,:,ismember(PCcellLines,cellLine)); 
                    objBounds                 = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured      = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);
                    [toxicDataset,geneEssentiality]   = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes  = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
                    filename = strcat('SGDpredictions/',cellLine,'_inPCpanel_',medium);
                    save(filename,'modelOpt','geneEssentiality','levyReport')
                case 'PC3'
                    load Models/cPC3.mat
                    load COREforLeviathan
                    refModel     = addDiffusionCORErxns(PC3);
                    refModel     = convertBoundsToFinite(refModel,1000);
                    obj          = {'CANCERBOUT'};
                    [~,metFluxes,metList]     = retrieveCommonExchangeFluxes(refModel,cellLine,2,obj,true,true,medium);
                    objBounds                 = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured      = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);
                    [toxicDataset,geneEssentiality]   = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes  = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
                    filename = strcat('SGDpredictions/',cellLine,'_inPCpanel_',medium);
                    save(filename,'modelOpt','geneEssentiality','levyReport')
                case '786Ov1'
                    load Models/786O_v1.mat
                    refModel = addDiffusionCORErxns(VHLmodel);
                    refModel = convertBoundsToFinite(refModel,1000);
                    RCCcellLines = unique(CLs(strncmp('Renal',Histology,5)));
                    obj   = {'CANCERBOUT'};
                    [~,allMetFluxes,metList] = retrieveCommonExchangeFluxes(refModel,RCCcellLines,2,obj,true,true,medium);
                    metFluxes    = allMetFluxes(:,:,ismember(RCCcellLines,cellLine)); 
                    objBounds    = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);                    
                    [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
                case '786Ov2'
                    load Models/786O_v2.mat
                    refModel = addDiffusionCORErxns(VHLmodel);
                    refModel = convertBoundsToFinite(refModel,1000);
                    RCCcellLines = unique(CLs(strncmp('Renal',Histology,5)));
                    obj   = {'CANCERBOUT'};
                    [~,allMetFluxes,metList] = retrieveCommonExchangeFluxes(refModel,RCCcellLines,2,obj,true,true,medium);
                    metFluxes    = allMetFluxes(:,:,ismember(RCCcellLines,cellLine)); 
                    objBounds    = [0,24*log(2)./mean(DT(ismember(CLs,cellLine)))]*550;
                    metMediumNotMeasured = setdiff(mediumMets,metList);
                    metMediumNotMeasuredBound = repmat([-1000 1000],numel(metMediumNotMeasured),1);
                    [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,...
                        [metList;metMediumNotMeasured],[metFluxes;metMediumNotMeasuredBound],...
                        objBounds,typical,parameters);     
                    [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold);
                    unviableGenes = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes); 
                    writeTxtFromCell(unviableGenes,strcat('SGDpredictions/unviableSGDin',refModelName,wFluxes,'_',cellLine,'_',medium));
                    drawPerturbationSummary(geneEssentiality,'sgd',strcat('Images/SGDin',refModelName,'_w',cellLine,'Fluxes_',medium))
            end
    end
end
%% APPENDIX: support functions
function modelCORE = addDiffusionCORErxns(model)
    diffusions.rxns = {'newDiffusion1';
    'newDiffusion2';
    'newDiffusion3';
    'newDiffusion4';
    'newDiffusion5';
    'newDiffusion6'};
    diffusions.equations = {'betaine[c] => betaine[s]';
        'citrate[c] => citrate[s]';
        'DHAP[c] => DHAP[s]';
        'GMP[c] => GMP[s]';
        'succinate[c] => succinate[s]';
        'thiamin[c] => thiamin[s]'};
    diffusions.subSystems = {'Fake reactions';
        'Fake reactions';
        'Fake reactions';
        'Fake reactions';
        'Fake reactions';
        'Fake reactions'};
    diffusions.lb   = zeros(6,1);
    diffusions.ub   = Inf*ones(6,1);
    if any(~ismember(diffusions.rxns,model.rxns))
        diffusionsToAdd = subsetstruct(diffusions,~ismember(diffusions.rxns,model.rxns));
        modelCORE       = addRxns(model,diffusionsToAdd,3); 
    else
        modelCORE = model;
    end
end
function s2        = subsetstruct(s1,is)
%SUBSETSTRUCT Extract subset of each field of a structure.
%    s2 = subsetstruct(s1,is)
%    s1 = full structure
%    is = n-element logical vector of field elements to extract, where n is
%         size of dimension to use for extraction (max of 2 dimensions, and
%         extracts first dimension if size of both dim are same--use 
%         function twice in succession in this case)
%       = if index vector, converted to n-element logical vector assuming
%         maximum dimension of first field as size n
%         (use logical vector if this is not the case)
%    s2 = sub structure

% Copyright (c) 1994-2013 by Michael G. Kay
% Matlog Version 15 07-Jan-2013 (http://www.ise.ncsu.edu/kay/matlog)

% Input Error Checking ****************************************************
if ~isstruct(s1) || length(s1) ~= 1
   error('"s1" must be a single structure.')
end
% End (Input Error Checking) **********************************************

fn = fieldnames(s1);
s2 = s1;
for i = 1:length(fn)
   fv = getfield(s1,fn{i});
   if i == 1
      if islogical(is)       %Logical
         n = length(is(:));
         idx = find(is); % Use index instead of logical so that can have
                         % multiple values the same (e.g., if multiple
                         % cites the same in idx, then s = uscity(idx)
                         % gives s the same size as idx; using logical
                         % only gives the unique values in idx)
      else                   %Index
         n = length(fv);
         %is = idx2is(is,n);
         idx = is;
      end
   end
   if size(fv,1) == n
      s2 = setfield(s2,fn{i},fv(idx,:)); % is -> idx
   elseif size(fv,2) == n
      s2 = setfield(s2,fn{i},fv(:,idx)); % is -> idx
   else
      error('Length of "is" does not match any field dimensions.')
   end
end
end
function FBSmets   = getOldFBSMets
    fid = fopen('FBScompositionInHMR.txt','r');
    t   = textscan(fid,'%s','Delimiter',';');
    fclose(fid);
    FBSmets = t{:};
end
function HAMmets   = getHAMMets
    fid = fopen('HAMcompositionInHMR.txt','r');
    t   = textscan(fid,'%s','Delimiter',';');
    fclose(fid);
    HAMmets = t{:};
end
function model     = fixMisannotatedGene(model)
    % Fix misannotated OGDH
    OGDH            = findEnsemblFromHugo('OGDH');   %Promiscous
    OGDHr           = getRxnsFromGenes(model,OGDH);
    wrongOGDHr      = OGDHr([1:7,9:10]);
    wrongOGDHrIndex = getIndexes(model,wrongOGDHr,'rxns');
    model.rxnGeneMat(wrongOGDHrIndex,getIndexes(model,OGDH,'genes')) = 0;
end
function [commonRxns,commonBounds,commonMets] = retrieveCommonExchangeFluxes(refModel,CL,relaxationTerm,obj,suppressWarnings,suppressOutput,medium)
        modelG = refModel;
        modelG = simplifyModel(modelG);
        modelG = addCompartmentalizedMetNames(modelG);
%         modelG = convertBoundsToInfinite(modelG,true);
        if ~suppressOutput
            fprintf('\t%s\t%s\n','Using typical constraints on uptake fluxes:','Jain et al. (2012)')
            fprintf('\t%s\t%d\n\t\t%s\n','N° of selected cell lines from the panel:',numel(CL))
            disp(CL)
        end
        % MATCH METABOLITES TO MODEL EXCHANGE FLUXES
        % Typical 2 is CORE profile from Jain et al., 2012. It 
        % includes 123 metabolites (note that the cells were cultured
        % in FBS medium).These will be searched against
        % the model metabolite list. If not found in the [s]
        % compartment, they will be searched in [c]. If still not
        % found, they will be ignored. If for some metabolites an 
        % exchange rxns could not be found, these will be added to
        % the model. Note: some extra exchange rxns  are added according
        % to evidence in the paper.
        if exist('COREforLeviathan.mat','file') ~= 2
            throw(MException('','The COREforLeviathan file containing the CORE flux profiles from Jain et al., 2012, is not in MATLAB search path'));
        else
            load COREforLeviathan.mat 
        end
        if ~all(ismember(CL,CLs))
            throw(MException('','Input cell line name(s) not found in the dataset from Jain et al., 2012.'));
        end
        modelDiff = addDiffusionCORErxns(modelG);
        if strcmp(medium,'FBS')
            FBSmets = getOldFBSMets;
            if ~suppressOutput
                fprintf('\t%s\n','Constraining the model to a Fetal Bovine Serum medium...')
            end
        elseif strcmp(medium,'HAM')
            FBSmets = getHAMMets;
            if ~suppressOutput
                fprintf('\t%s\n','Constraining the model to a HAM medium...')
            end
        end
        if ~all(ismember(FBSmets,modelDiff.metNamesC))
            FBSmets = FBSmets(ismember(FBSmets,modelDiff.metNamesC));
            if ~suppressWarnings
                fprintf('\t%s\n','Warning: not all mets in the medium are included in the model. These will be ignored');
            end
        end
        [~,metUptRxnsFBS,modelUptRxnsFBS]   = getUptakeReactionsMatchedToMetList(modelDiff,FBSmets,suppressWarnings,suppressOutput);
        modelFBS                            = constrainUptRxnsToFluxes(modelUptRxnsFBS,metUptRxnsFBS);
        % CONSTRAIN EXCHANGE FLUXES
        % Typical 2 constrains the model to CORE fluxes. Most likely,
        % these constraints will render the model unfeasible. By default,
        % the script removes all unfeasible constraints according to
        % the minimum global amount of individual fluxes required to 
        % fit the CORE flux profile. 
        % Note: the script starts evaluating the first cell line.
        nCLs       = numel(CL);
        selectedCL = CL(1);      
        if ~suppressOutput
            fprintf('\t%s\n','Constraining the model to a Jain et al. (2012) measured fluxes...')
        end  
        if ~all(ismember(CORE_mets_adj,modelFBS.metNamesC))
            CORE_mets_adj = CORE_mets_adj(ismember(CORE_mets_adj,modelFBS.metNamesC));
            if ~suppressWarnings
                fprintf('\t%s\n','Warning: not all mets in CORE are included in the model. These will be ignored');
            end
        end
        [~,metUptRxnsCORE,modelUptRxnsCORE] = getUptakeReactionsMatchedToMetList(modelFBS,CORE_mets_adj,suppressWarnings,suppressOutput);
        if ~suppressOutput
            fprintf('\t%s\n','Diagnosing whether imposed constraints are feasible...')
        end
        [~,commonRxns] = setCOREboundsForCL(modelUptRxnsCORE,selectedCL,metUptRxnsCORE,relaxationTerm,true,suppressWarnings,suppressOutput);
        [~,commonRxns] = setCOREboundsForCL(modelUptRxnsCORE,selectedCL,commonRxns,relaxationTerm,true,suppressWarnings,suppressOutput);
        if numel(CL) > 1
            for i = 2:nCLs
                selectedCL = CL(i);                
                [~,feasRxns]   = setCOREboundsForCL(modelUptRxnsCORE,selectedCL,metUptRxnsCORE,relaxationTerm,true,suppressWarnings,suppressOutput);
                [~,feasRxns]   = setCOREboundsForCL(modelUptRxnsCORE,selectedCL,feasRxns,relaxationTerm,true,suppressWarnings,suppressOutput);
                commonRxns = intersect(commonRxns,feasRxns);
            end
        end              
        % CHECK FOR VIABILITY
        % The model must be viable prior to any gene perturbation.
        % Since no constraints is a finite number, an artificial
        % upper bound for the objective function is given (i.e. 1).
        for i = 1:nCLs      
            selectedCL = CL(i);
            [modelCORE,commonRxns]  = setCOREboundsForCL(modelUptRxnsCORE,selectedCL,commonRxns,relaxationTerm,false,suppressWarnings,suppressOutput);
            commonMets          = getMetListMatchedToUptakeReactions(modelG,commonRxns);
            commonBounds(:,:,i) = [modelCORE.lb(getIndexes(modelCORE,commonRxns,'rxns')),modelCORE.ub(getIndexes(modelCORE,commonRxns,'rxns'))];
            CLind      = find(strcmp(selectedCL,CLs)); 
            bounds     = [0,24*log(2)./mean(DT(CLind(1)))]*550;
            [optim(i),modelOpt(i)] = optimizeObjective(modelCORE,obj,bounds,1,suppressOutput); %#ok<AGROW>
            if optim(i).stat ~= 1 || optim(i).x(getIndexes(modelCORE,obj,'rxns')) < eps
                throw(MException('','The model could not be optimized with implemented objective and CORE constraints'));
            end
            modelOpt(i).description = selectedCL; %#ok<AGROW>
        end            
        % Congrats! A condition-specific model(s) has been created. 
        % Jump to part B (apply gene perturbation).   
end
function [modelCORE,retainedRxns]  = setCOREboundsForCL(model,CL,metUptRxnsCORE,relaxationTerm,fixRigidBounds,suppressWarnings,suppressOutput)   
        load COREforLeviathan.mat        
        if isempty(metUptRxnsCORE)
            [~,metUptRxnsCORE] = getUptakeReactionsMatchedToMetList(model,CORE_mets_adj,suppressWarnings,suppressOutput);
        end
        CLinds = find(strcmp(CL,CLs));        
        % Detect bounds in CORE matrix and fix them. The first step
        % reassigns the metabolite list that indexes the matched uptake
        % reactions to the matrix, given that some metabolites might have
        % switched compartment when uptake reactions were searched or might
        % have been deleted because they were dead end in the model
        metListCORE    = getMetListMatchedToUptakeReactions(model,metUptRxnsCORE);
        metListCORE_nc = cellfun(@(a) a(1:end-3),metListCORE,'uni',false);
        CORE_mets_nc   = cellfun(@(a) a(1:end-3),CORE_mets_adj,'uni',false);
        [~,w]          = ismember(metListCORE_nc,CORE_mets_nc);
        fluxValues     = -CORE_adj(w,CLinds);
        meanFlux       = mean(fluxValues,2);
        stdFlux        = std(fluxValues,0,2);
            % Some metabolites were not measured in duplicate. The standard
            % deviation is estimated as the average standard deviation
            % among all the remaining metabolites
            stdFlux(stdFlux==0) = mean(stdFlux(stdFlux>0));
        %Equality constraints are replaced by relaxed upper and lower
        %bounds that are dependant on the standard deviation
        lbFlux      = meanFlux - relaxationTerm*stdFlux;  
        ubFlux      = meanFlux + relaxationTerm*stdFlux;
        modelCORE   = setParam(model,'lb',metUptRxnsCORE,lbFlux);
        modelCORE   = setParam(modelCORE,'ub',metUptRxnsCORE,ubFlux);
        if fixRigidBounds
            [modelCORE,retainedRxns] = fixUnfeasibleConstraints(modelCORE,metUptRxnsCORE,suppressWarnings,suppressOutput);
        else
            retainedRxns   = metUptRxnsCORE;
        end
end
function modelFin   = convertBoundsToFinite(model,arbitraryBound)
% convertBoundsToFinite
%
% This function converts all infinite bounds to +/-arbitraryBound. 
% If no arbitraryBound is passed, it will be assigned the largest absolute
% bound in the model.
%
% INPUT
% model:    a genome-scale metabolic model structure.
% arbitraryBound: a double containing the arbitrary value to assign in
%                 place of the infinite bounds.
% OUTPUT
% modelFin: a genome-scale metabolic model structure with bounds converted
%           to a finite number.
%
% Usage:  modelFin   = convertBoundsToFinite(model,arbitraryBound)
%
% 2013-10-23 Francesco Gatto
        if nargin < 2
            tmpModel = model;
            tmpModel.ub(tmpModel.ub==+Inf)=0;            
            tmpModel.lb(tmpModel.lb==-Inf)=0;
            arbitraryBound = max(max(tmpModel.ub),max(abs(tmpModel.lb)));
        end        
        model.ub(model.ub==+Inf)=+arbitraryBound;
        model.lb(model.lb==-Inf)=-arbitraryBound;        
        modelFin = model;
end
function [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,metList,metFluxes,objBounds,typical,parameters,tasks)
% function assessModelViability
%   This is the main function of the Leviathan toolbox. The function takes
%   a genome-scale model and checks for its viability upon gene
%   perturbations in a condition specified by the user. The reference model 
%   is constrained to a certain lifestyle determined by user-defined flux 
%   profiles of extracellular metabolites.
%   The script consists of 3 parts:
%   A) Create a condition-specific model;
%   B) Apply a panel of gene perturbations;
%   C) Validate lethal perturbations.
%
% INPUT
%   refModel:   a genome-scale model that will be tailored upon the
%               chosen test to be either constrained to get certain fluxes
%               or to selectively express certain reactions.
%   obj:    a cell array with the reaction IDs that are supposed to be
%           optimized (def: max) by the model to be considered viable.
%   perturbation: a string indicating which perturbation is to apply,
%                 either 'sgd' (single-gene deletion), 'dgd'
%                 (double-gene deletion), or 'am' (antimetabolite).
%   metList:  either a cell array of metabolites or the filename of a plain
%       text file containing metabolite names in each row. Note that
%       metabolites can be either compartmentalized or not (e.g. glucose OR
%       glucose[c]).       
%   metFluxes: a nx1 or nx2 numerical vector array of exchange fluxes for
%           each metabolite in metList. If nx2 is provided, then
%           the first column is assumed to be a vector of lower bounds,
%           and the second as upper bounds (default: -x to +x, where x is
%           the largest absolute value in the model.lb | model.ub).
%   objBounds: only if typical is not selected. A nx1 or nx2 numerical 
%              vector array for the upper and lower bound used for the 
%              objective flux(es) specified in obj. If nx1 is provided,
%              an equality constraint is assumed (default: [0,1]).
%   typical: a logical that defines whether a typical medium 
%            should be used to constrain lifestyle (Fetal Bovine Serum)
%   parameters: a structure containing the following fields:
%       - suppressWarnings: a logical if warnings on screen should be
%                     suppressed (default: false)
%       - suppressOutput: a logical if output on screen should be
%                     suppressed (default: true)
%       - minmax:   a numerical, where 1 is maximization and -1 is
%                   minimization (default: 1).
%       - useControl: a logical if only perturbations that are not
%                     unviable to the reference model should be scanned
%                     in the condition-specific model (default: false).
%       - relaxationTerm: a numerical indicating the percent of
%                         relaxation on the bounds imposed via
%                         metFluxes (only if metFluxes is nx2, def: 
%                         10*eps).
%       - fixUnfeasibleFluxes: a logical indicating whether the minimum amount of
%                         constraints should be removed in order to get a feasible
%                         solution (i.e. steady-state flux-balanced
%                         distribution that optimizes the objective)
%                         (default: true). Warning: if false, make sure that
%                         the bounds in "metFluxes" are feasible before
%                         running the script.
%       - genePerturbationMethod: a string indicating whether 'fba' or
%                             'moma' should be used when applying 'sgd'
%                             or 'dgd' (def: 'fba').
%       - singleGenesNotToPerturb: a nx1 cell array of string containing a list
%                            of genes that should not be perturbed when
%                            applying either 'sgd' or 'dgd' (def:
%                            empty).
%       - doubleGenesNotToPerturb: a nx2 cell array of string containing a list
%                            of gene pairs that should not be perturbed when
%                            applying 'dgd' (def: empty).
%       - excludeSingleGenesNotToPerturbFromDgd: a logical indicating
%                         whether the list of genes to pair to test DGD
%                         should be obtained by all the genes in the
%                         model except those indicated in
%                         "singleGenesNotToPerturb" (def: false). Only 
%                         if "perturbation" is 'dgd'. Note: it can be used
%                         even if "singleGenesNotToPerturb" is empty, in
%                         this case the script will perform all possible
%                         double gene perturbations.
%       - convertBoundsToInfinite: a logical if arbitrarily large bounds
%                                  should be replaced by +/-Inf. This is
%                                  not recommended if few constraints are
%                                  imposed (default: false).%  
%       - saveWTmodel: a logical or a string if the condition-specific model 
%                      that can be optimized before applying any gene knock-out 
%                      should be saved. If string, this defines the
%                      filename (opt, default: false).
%       - saveSGD: a logical or a string if the SGD results
%                      should be saved. If string, this defines the
%                      filename (opt, default: false).
%   tasks: a task structure obtained with parseTaskList or a logical 
%          indicating whether a list of common metabolic tasks that 
%          should be fulfilled by the model (def: false). The tasks are
%          checked upon each individual perturbation in the
%          condition-specific model.
% OUTPUT
%   levyReport: a structure containing the following fields
%       - addedRxns: a cell array of rxn IDs corresponding to those rxns
%                    that were added to the reference model in order to
%                    accomodate the uptake of the metabolites in the input.
%       - feasMets: a cell array of metabolites that have been successfully
%                   among those provided in metList (only if
%                   fixUnfeasibleFluxes is true).
%       - feasRxns: a cell array of rxn IDs that were constrained
%                   matching feasMets (only if fixUnfeasibleFluxes is true).
%       - unviableInControl: a cell array of single perturbations unviable
%                            in a control model (if available).
%       - unviableDGDinControl: a cell array of double gene deletions unviable
%                            in a control model (if available).
%       - taskReport: a structure containing the following fields:
%           - taskFeasibility: a n x o binary matrix, where each entry
%                              is either 1 or 0 if n-th task given the o-th
%                              perturbation could be accomplished or not.
%           - taskNames: a cell array of strings that indicate the
%                        metabolic tasks checked for feasibility for each
%                        perturbation. It indexes the n dimension in the 
%                        taskFeasibility matrix.
%           - SGD: only if perturbation is 'sgd'. A cell array of strings
%                  that indicate the lethal perturbations tested for tasks.
%                  It indexes the o dimension in the taskFeasibility
%                  matrix.
%           - DGD: only if perturbation is 'dgd'. A cell array of strings
%                  that indicate double lethal perturbations tested for tasks.
%                  It indexes the o dimension in the taskFeasibility
%                  matrix.
%
% Usage: [modelOpt,geneEssentiality,levyReport] = assessModelViability(refModel,obj,perturbation,metList,metFluxes,objBounds,typical,parameters,tasks)
%
% Francesco Gatto, 2014-08-26

%% ARGUMENT CHECK
% Check inputs
if nargin < 9  || isempty(tasks)
    checkTasks = false;
else
    if ~isstruct(tasks) 
        if ~islogical(tasks)
            throw(MException('','If tasks are to be checked, feed either a task-containing structure or a logical to supply a list of common metabolic tasks'));
        else
            if tasks
                load TasksForLeviathan.mat
                checkTasks = true;
            else
                checkTasks = false;
            end
        end
    else
        checkTasks = true;
    end
end
if nargin < 8 || isempty(parameters)
    suppressWarnings    = false;
    suppressOutput      = true;
    minmax              = 1;
    useControl          = false;
    fixUnfeasibleFluxes = true;
    relaxationTerm      = [];
    genePerturbationMethod  = 'fba';
    singleGenesNotToPerturb = [];
    doubleGenesNotToPerturb = [];
    excludeSingleGenesNotToPerturbFromDgd = false;
    convertBoundsToInfinite = false;
    saveWTmodel = false;
    saveSGD     = false;
else
    [suppressWarnings,suppressOutput,minmax,useControl,relaxationTerm,fixUnfeasibleFluxes,...
    genePerturbationMethod,singleGenesNotToPerturb,doubleGenesNotToPerturb,...
    excludeSingleGenesNotToPerturbFromDgd,convertBoundsToInfinite,saveWTmodel,...
    saveWTmodelFilename,saveSGD,saveSGDfilename] = checkParameters(parameters);
end
if nargin < 7 || isempty(typical)
    typical = false;
else
    if ~islogical(typical)
        throw(MException('','Typical must be a logical'));
    end
end
if nargin < 6 || isempty(objBounds)
    objBounds = [0,1];
else
    if ~isnumeric(objBounds)
        throw(MException('','objBounds must be a 1x1 or 1x2 vector'));
    end
end      
if nargin < 5
    throw(MException('','A list of flux values to constrain the metabolites in "metList" must be given if typical is false'));
else
    if ~typical
        if isempty(metFluxes)
            throw(MException('','A list of flux values to constrain the metabolites in "metList" must be given if typical is false'));
        else
            if size(metFluxes,2) ~= 1 && size(metFluxes,2) ~= 2        
                throw(MException('','The list of flux values to constrain the metabolites in "metList" must be a nx1 or nx2 matrix'));
            end
            if size(metFluxes,1) ~= numel(metList)
                throw(MException('','The number of flux values to constrain the metabolites in "metList" does not equal the number of metabolites in in "metList"'));
            end
        end
    else
        if ~isempty(metFluxes)
            if size(metFluxes,2) ~= 1 && size(metFluxes,2) ~= 2        
                throw(MException('','The list of flux values to constrain the metabolites in "metList" must be a nx1 or nx2 matrix'));
            end
            if size(metFluxes,1) ~= numel(metList)
                throw(MException('','The number of flux values to constrain the metabolites in "metList" does not equal the number of metabolites in in "metList"'));
            end
        end
    end
end   
if nargin < 4
    throw(MException('','A list of metabolites whose exchange flux must be constrained must be given if typical is false'));
else
    if ~typical
        if isempty(metList)
            throw(MException('','A list of metabolites whose exchange flux must be constrained must be given if typical is false'));
        else
            if ischar(metList)
                fid  = fopen(metList,'r');
                if fid < 0
                    throw(MException('','The file containing the metabolite whose exchange flux should be constrained could not be found or is not well formatted (one metabolite per row)'));
                end
                scan = textscan(fid,'%s','Delimiter','\n');
                fclose(fid);
                metList = scan{1};
            end
            if ~iscell(metList)
                throw(MException('','The list of metabolites whose exchange flux must be constrained must be a cell array of strings'));
            else
                if numel(metList) ~= numel(unique(metList))
                    throw(MException('','The list of metabolites whose exchange flux must be constrained is not unique'));
                end
            end
            if all(cellfun(@(a) strcmp(a(end),']'),metList)) && all(cellfun(@(a) strcmp(a(end-2),'['),metList))
                isCompart = true;
            elseif ~any(cellfun(@(a) strcmp(a(end),']'),metList)) || ~any(cellfun(@(a) strcmp(a(end-2),'['),metList))
                isCompart = false;
            else
                throw(MException('','metList is expected to be either all compartmentalized (e.g. glucose[s]) or any of them'));
            end   
        end
    else
        if ~isempty(metList)
            if ischar(metList)
                fid  = fopen(metList,'r');
                if fid < 0
                    throw(MException('','The file containing the metabolite whose exchange flux should be constrained could not be found or is not well formatted (one metabolite per row)'));
                end
                scan = textscan(fid,'%s','Delimiter','\n');
                fclose(fid);
                metList = scan{1};
            end
            if ~iscell(metList)
                throw(MException('','The list of metabolites whose exchange flux must be constrained must be a cell array of strings'));
            else
                if numel(metList) ~= numel(unique(metList))
                    throw(MException('','The list of metabolites whose exchange flux must be constrained is not unique'));
                end
            end
            if all(cellfun(@(a) strcmp(a(end),']'),metList)) && all(cellfun(@(a) strcmp(a(end-2),'['),metList))
                isCompart = true;
            elseif ~any(cellfun(@(a) strcmp(a(end),']'),metList)) || ~any(cellfun(@(a) strcmp(a(end-2),'['),metList))
                isCompart = false;
            else
                throw(MException('','metList is expected to be either all compartmentalized (e.g. glucose[s]) or any of them'));
            end     
        end
    end
end
if ~suppressOutput
    fprintf('%s\n%s\n%s\n','%%%%','Initialization - Started','%%%%')
end
% Check objective
if ~iscell(obj) 
    if ~ischar(obj)        
        throw(MException('','Objective function must be either a string or a cell array of strings.'));
    else
        obj = {obj};
    end
else
    if ~ismember(obj,refModel.rxns)
        throw(MException('','Objective function not found in reference model.'));
    end
end
% Check perturbation
if ~ischar(perturbation) 
    throw(MException('','Perturbation must be a string, either "sgd","dgd", or "am"'));
else
    if ~strcmp(perturbation,{'sgd','dgd','am'})
        throw(MException('','Perturbation must be a string, either "sgd","dgd", or "am"'));
    end
end
% Check model
if ~isstruct(refModel)
        throw(MException('','Reference model must be a valid RAVEN formatted structure'));
else
    if ~isfield(refModel,'S')
        throw(MException('','Model structure must contain a stoichiometric matrix field named "S"'));
    else
        [nRefMets,nRefRxns] = size(refModel.S);
    end
    if ~isfield(refModel,'rxnGeneMat')
        throw(MException('','Model structure must contain a rxn-gene association matrix field named "rxnGeneMat"'));
    else
        [nRefRxns_red,nRefGenes] = size(refModel.rxnGeneMat);
        if numel(nRefRxns_red) ~= numel(nRefRxns)
            throw(MException('','N° of rxns in the rxn-gene association matrix is different than in the stoichiometric matrix.'));
        end
    end
end
% Report inputs
if ~suppressOutput
    fprintf('%s\n%s\n','Starting Leviathan Toolbox','Checking inputs...')
    if isfield(refModel,'description')
        fprintf('\t%s\t%s\n','Reference model loaded: ',refModel.description);
    else
        fprintf('\t%s\n','Reference model loaded - No description available');
    end
    fprintf('\t\t%s%d\n\t\t%s%d\n\t\t%s%d\n','N° rxns in the reference model: ',nRefRxns,...
        'N° mets in the reference model: ',nRefMets,'N° genes in the reference model: ',nRefGenes);
    fprintf('\t%s\t%s\n','Selected perturbation: ',perturbation);
    objEqn = constructEquations(refModel,obj);
    fprintf('\t%s\t%s\n\t\t%s\n','Selected objective: ',obj{:},objEqn{:});
    fprintf('\t%s\n\t\t%s\t%d\n\t\t%s\t%d\n\t\t%s\t%d\n\t\t%s\t%s\n\t\t%s\t%d\n\t\t%s\t%d\n\t\t%s\t%d\n','Selected parameters: ','suppressWarnings: ',suppressWarnings,...
        'useControl: ',useControl,'fixUnfeasibleFluxes: ',fixUnfeasibleFluxes,'genePerturbationMethod: ',genePerturbationMethod,...
        'singleGenesNotToPerturb: ',numel(singleGenesNotToPerturb),'doubleGenesNotToPerturb: ',numel(doubleGenesNotToPerturb)/2,...
        'excludeSingleGenesNotToPerturbFromDgd: ',excludeSingleGenesNotToPerturbFromDgd);
    fprintf('\t%s\t%d\n','Check Tasks? ',checkTasks);
end     
%% PREPROCESSING
% Preprocess refModel a bit. First simplify it (such that exchange rxns are
% in the form ' <=> A'), then add the field of compartmentalized
% metabolites, next convert all finite bounds to infinite.
if ~suppressOutput
    fprintf('%s\n','Preprocessing reference model...')
end
modelG = refModel;
modelG = simplifyModel(modelG);
modelG = addCompartmentalizedMetNames(modelG);
if convertBoundsToInfinite
    modelG       = convertBoundsToInfinite(modelG,true);
    largestBound = Inf;
else
    modelG       = convertBoundsToFinite(modelG);
    largestBound = max(max(modelG.ub),max(abs(modelG.lb)));
end
if ~suppressOutput
    fprintf('\t%s\n','Simplified model:')
    fprintf('\t\t%s%d\n\t\t%s%d\n\t\t%s%d\n','N° rxns in the simplified model: ',size(modelG.S,2),...
        'N° mets in the simplified model: ',size(modelG.S,1),'N° genes in the simplified model: ',size(modelG.rxnGeneMat,2));
    fprintf('\t%s\t%d\t%s%s%s\n','Added compartmentalized metabolite name field:',numel(modelG.metNamesC),' (e.g. ',modelG.metNamesC{1},')')
    fprintf('\t%s\t%d\n','Artificial maximum and minimim bounds were converted to:',largestBound)
    fprintf('%s\n%s\n%s\n\n','%%%%','Initialization - Completed','%%%%')
end
%% A) CREATE A CONDITION-SPECIFIC MODEL
if ~suppressOutput
    fprintf('%s\n%s\n%s\n','%%%%','Part A - Started','%%%%')
    fprintf('%s\n','Creating a condition-specific model...')
end
if typical
    if ~suppressOutput
        fprintf('\t%s\n','Constraining exchange fluxes to typical medium')
    end
    FBSmets              = getFBScomposition;
    [~,FBSUptRxns,modelG]= getUptakeReactionsMatchedToMetList(modelG,FBSmets,suppressWarnings,suppressOutput);
end
if ~suppressOutput
    fprintf('\t%s\n','Constraining exchange fluxes to user-defined inputs')
end
[acceptedMetList,acceptedUptRxns,modelUptRxns] = getUptakeReactionsMatchedToMetList(modelG,metList,suppressWarnings,suppressOutput);
if ~isempty(metList)
    if isCompart
        acceptedMetFluxes = metFluxes(ismember(cellfun(@(a) a(1:end-3),metList,'uni',false),cellfun(@(a) a(1:end-3),acceptedMetList,'uni',false)),:);
    else    
        acceptedMetFluxes = metFluxes(ismember(metList,cellfun(@(a) a(1:end-3),acceptedMetList,'uni',false)),:);
    end
else
    acceptedMetFluxes = [];
end
if typical
    modelConstrained  = constrainUptRxnsToFluxes(modelUptRxns,[acceptedUptRxns;setdiff(FBSUptRxns,acceptedUptRxns)],[acceptedMetFluxes;repmat([-largestBound largestBound],numel(setdiff(FBSUptRxns,acceptedUptRxns)),1)],relaxationTerm);
    constraints       = dataset(getMetListMatchedToUptakeReactions(modelUptRxns,[acceptedUptRxns;setdiff(FBSUptRxns,acceptedUptRxns)]),[acceptedMetFluxes;repmat([-largestBound largestBound],numel(setdiff(FBSUptRxns,acceptedUptRxns)),1)],'VarNames',{'Metabolites';'Bounds_on_Uptake_rxn'});
else
    modelConstrained  = constrainUptRxnsToFluxes(modelUptRxns,acceptedUptRxns,acceptedMetFluxes,relaxationTerm);
    constraints       = dataset(getMetListMatchedToUptakeReactions(modelUptRxns,acceptedUptRxns),acceptedMetFluxes,'VarNames',{'Metabolites';'Bounds_on_Uptake_rxn'});
end
if ~isempty(metList)
    if fixUnfeasibleFluxes
        if ~suppressOutput
           fprintf('\t\t%s\n','Diagnosing whether imposed constraints are feasible...')
        end
        [modelConstrained,feasRxns] = fixUnfeasibleConstraints(modelConstrained,acceptedUptRxns,suppressWarnings,suppressOutput);
        levyReport.feasMets  = getMetListMatchedToUptakeReactions(modelConstrained,feasRxns);
        levyReport.feasRxns  = feasRxns;                
    end
end    
[optim,modelOpt] = optimizeObjective(modelConstrained,obj,objBounds,minmax,[],suppressOutput);
if optim.stat ~= 1 || optim.x(getIndexes(modelOpt,obj,'rxns')) < eps
    throw(MException('','The model could not be optimized with implemented objective and user-defined fluxes'));
end
if saveWTmodel
    save(saveWTmodelFilename,'modelOpt')
end
levyReport.addedRxns   = setdiff(modelOpt.rxns,refModel.rxns);
levyReport.constraints = constraints;
if ~suppressOutput
    fprintf('%s\n%s\n%s\n\n','%%%%','Part A - Completed','%%%%')
end
%% B) CHECK FOR VIABILITY OF THE MODEL UPON PERTURBATION
% This part revolves around the function applyGenePerturbation
if ~suppressOutput
    fprintf('%s\n%s\n%s\n','%%%%','Part B - Started','%%%%')
    fprintf('%s\n','Applying perturbation...')
end
switch perturbation
    case 'sgd'
        settings.genePerturbationMethod  = genePerturbationMethod;
        settings.singleGenesToExclude    = singleGenesNotToPerturb;
        settings.suppressOutput          = suppressOutput;
        if saveSGD
            settings.saveSGDfilename     = saveSGDfilename;
        end
        if useControl
            if ~suppressOutput
                fprintf('\t%s\n','In control model:')
            end
            [~,modelGopt]    = optimizeObjective(modelG,obj,[0,1]);
            geneEss_control  = applyGenePerturbation(modelGopt,perturbation,settings);
            unviableInControl= union(geneEss_control.lethalGenes,geneEss_control.toxicGenes);
            singleGenesToExclude = [singleGenesNotToPerturb;unviableInControl];
            settings.singleGenesToExclude = singleGenesToExclude;
        end
        if ~suppressOutput
            fprintf('\t%s\n','In condition-specific model:')
        end
        geneEssentiality = applyGenePerturbation(modelOpt,perturbation,settings);
        if useControl
            levyReport.unviableInControl = unviableInControl; 
        end
    case 'dgd'
        % It is usually smart not to perturb a priori SGD found essential in
        % control and case by applying this pipeline for "sgd"
        settings.genePerturbationMethod  = genePerturbationMethod;
        settings.singleGenesToExclude = singleGenesNotToPerturb;
        settings.doubleGenesToExclude = doubleGenesNotToPerturb;
        settings.excludeSingleGenesNotToPerturbFromDgd = excludeSingleGenesNotToPerturbFromDgd;
        settings.suppressOutput          = suppressOutput;
        if useControl
            if ~suppressOutput
                fprintf('\t%s\n','In control model:')
            end
            [~,modelGopt]          = optimizeObjective(modelG,obj,[0,1]);
            geneEss_control        = applyGenePerturbation(modelGopt,perturbation,settings);
            singleUnviableInControl= geneEss_control.lethalSingleGenes;
            singleGenesToExclude   = [singleGenesNotToPerturb;singleUnviableInControl];
            doubleUnviableInControl= [geneEss_control.lethalDoubleGenes;geneEss_control.toxicDoubleGenes];
            doubleGenesToExclude   = [doubleGenesNotToPerturb;doubleUnviableInControl];
            settings.singleGenesToExclude = singleGenesToExclude;
            settings.doubleGenesToExclude = doubleGenesToExclude;
        end
        if ~suppressOutput
            fprintf('\t%s\n','In condition-specific model:')
        end
        geneEssentiality = applyGenePerturbation(modelOpt,perturbation,settings);
        if useControl
            levyReport.unviableDGDInControl = doubleUnviableInControl; 
        end
    case 'am'
        return
end
if ~suppressOutput
    fprintf('%s\n%s\n%s\n','%%%%','Part B - Completed','%%%%')
end
%% C) VALIDATE CONDITION-SPECIFIC ESSENTIALITY
% This part revolves around the function getTaskReport applied to proper
% perturbations via applyGeneKO
if checkTasks
    switch perturbation
        case 'sgd'
                lethalPerturbations  = union(geneEssentiality.lethalGenes,geneEssentiality.toxicGenes);
                nLethalPerturbations = numel(lethalPerturbations);
                nTasks               = numel(tasks);
                taskFeasibility      = zeros(nTasks+1,nLethalPerturbations);
                for i = 1:nLethalPerturbations
                    modelKO = applyGeneKO(modelOpt,lethalPerturbations(i));
                    [taskFeasibility(:,i),taskNames] = getTaskReport(modelKO,tasks);
                end
                taskReport.taskNames = taskNames;
                taskReport.SGD       = lethalPerturbations;
                taskReport.taskFeasibility = taskFeasibility;
        case 'dgd'
            % Dgd differ from sgd in that the output from part B) is a nx2
            % cell array of gene pairs thare are unviable only in the
            % condition-specific model.
                lethalPerturbations  = unviableDGDonlyInCase;
                nLethalPerturbations = size(lethalPerturbations,1);
                nTasks               = numel(tasks);
                taskFeasibility      = zeros(nTasks,nLethalPerturbations);
                for i = 1:nLethalPerturbations
                    modelKO = applyGeneKO(modelOpt,lethalPerturbations(i,:));
                    [taskFeasibility(:,i),taskNames] = getTaskReport(modelKO,tasks);
                end
                taskReport.taskNames = taskNames;
                taskReport.DGD       = lethalPerturbations;
                taskReport.taskFeasibility = taskFeasibility;
        case 'am'
            return
    end
    levyReport.taskReport = taskReport;
end
end
function [suppressWarnings,suppressOutput,minmax,useControl,relaxationTerm,fixUnfeasibleFluxes,...
    genePerturbationMethod,singleGenesNotToPerturb,doubleGenesNotToPerturb,...
    excludeSingleGenesNotToPerturbFromDgd,convertBoundsToInfinite,saveWTmodel,saveWTmodelFilename,...
    saveSGD,saveSGDfilename] = checkParameters(parameters)
    if isfield(parameters,'suppressWarnings')
        suppressWarnings = parameters.suppressWarnings;
        if ~islogical(suppressWarnings)
            disp('Warning: suppressWarnings is not a logical. Default opt (false) is set\n')
            suppressWarnings = false;
        end
    else
        suppressWarnings = false;
    end
    if isfield(parameters,'suppressOutput')
        suppressOutput = parameters.suppressOutput;
        if ~islogical(suppressOutput)
            if  ~suppressWarnings
                disp('Warning: suppressOutput is not a logical. Default opt (true) is set\n')
                suppressOutput = true;
            end            
        end
    else
        suppressOutput = true;
    end
    if isfield(parameters,'minmax')
        minmax = parameters.minmax;
        if ~isnumeric(minmax)
            if  ~suppressWarnings
                disp('Warning: minmax is not a numeric. Default opt (1: max) is set\n')
            end
            minmax = 1;
        end
    else
        minmax = 1;
    end
    if isfield(parameters,'useControl')
        useControl = parameters.useControl;
        if ~islogical(useControl)
            if  ~suppressWarnings
                disp('Warning: useControl is not a logical. Default opt (true) is set\n')
            end
            useControl = true;
        end
    else
        useControl = false;
    end
    if isfield(parameters,'relaxationTerm')
        relaxationTerm = parameters.relaxationTerm;
        if ~isnumeric(relaxationTerm) && ~isempty(relaxationTerm);
            if ~suppressWarnings
                disp('Warning: relaxationTerm is not numeric. Default opt (empty) is set\n')
            end
            relaxationTerm = [];
        end
    else
        relaxationTerm = [];
    end
    if isfield(parameters,'fixUnfeasibleFluxes')
        fixUnfeasibleFluxes = parameters.fixUnfeasibleFluxes;
        if ~islogical(fixUnfeasibleFluxes)
            if  ~suppressWarnings
                disp('Warning: fixUnfeasibleFluxes is not a logical. Default opt (true) is set\n')
            end
            fixUnfeasibleFluxes = true;
        end
    else
        fixUnfeasibleFluxes = true;
    end
    if isfield(parameters,'genePerturbationMethod')
        genePerturbationMethod = parameters.genePerturbationMethod;
        if ~ischar(genePerturbationMethod)
            if  ~suppressWarnings
                disp('Warning: genePerturbationMethod is not a string. Default opt ("fba") is set\n')
            end
            genePerturbationMethod = 'fba';
        end
    else
        genePerturbationMethod = 'fba';
    end
    if isfield(parameters,'singleGenesNotToPerturb')
        singleGenesNotToPerturb = parameters.singleGenesNotToPerturb;
        if ~iscell(singleGenesNotToPerturb)
            if  ~suppressWarnings
                disp('Warning: singleGenesNotToPerturb is not a cell array. Default opt ("empty") is set\n')
            end
            singleGenesNotToPerturb = [];
        end
    else
        singleGenesNotToPerturb = [];
    end
    if isfield(parameters,'doubleGenesNotToPerturb')
        doubleGenesNotToPerturb = parameters.doubleGenesNotToPerturb;
        if ~iscell(doubleGenesNotToPerturb)
            if  ~suppressWarnings
                disp('Warning: doubleGenesNotToPerturb is not a cell array. Default opt ("empty") is set\n')
            end
            doubleGenesNotToPerturb = [];
        end
    else
        doubleGenesNotToPerturb = [];
    end
    if isfield(parameters,'excludeSingleGenesNotToPerturbFromDgd')
        excludeSingleGenesNotToPerturbFromDgd = parameters.excludeSingleGenesNotToPerturbFromDgd;
        if ~islogical(excludeSingleGenesNotToPerturbFromDgd)
            if  ~suppressWarnings
                disp('Warning: excludeSingleGenesNotToPerturbFromDgd is not a logical. Default opt (false) is set\n')
            end
            excludeSingleGenesNotToPerturbFromDgd = false;
        end
    else
        excludeSingleGenesNotToPerturbFromDgd = false;
    end
    if isfield(parameters,'convertBoundsToInfinite')
        convertBoundsToInfinite = parameters.convertBoundsToInfinite;
        if ~islogical(convertBoundsToInfinite)
            if  ~suppressWarnings
                disp('Warning: convertBoundsToInfinite is not a logical. Default opt (false) is set\n')
            end
            convertBoundsToInfinite = false;
        end
    else
        convertBoundsToInfinite = false;
    end
    if isfield(parameters,'saveWTmodel')
        saveWTmodel = parameters.saveWTmodel;
        if ~islogical(saveWTmodel)
            if ~ischar(saveWTmodel)
                if  ~suppressWarnings
                    disp('Warning: saveWTmodel is not a logical. Default opt (false) is set\n')
                end
                saveWTmodel = false;
            else
                saveWTmodelFilename = saveWTmodel;
                saveWTmodel = true;
            end
        else
            saveWTmodelFilename = 'modelWT.mat';
        end
    else
        saveWTmodel = false;
    end
    if isfield(parameters,'saveSGD')
        saveSGD = parameters.saveSGD;
        if ~islogical(saveSGD)
            if ~ischar(saveSGD)
                if  ~suppressWarnings
                    disp('Warning: saveSGD is not a logical. Default opt (false) is set\n')
                end
                saveSGD = false;
            else
                saveSGDfilename = saveSGD;
                saveSGD = true;
            end
        else
            saveSGDfilename = 'modelSGD.txt';
        end
    else
        saveSGD = false;
    end
end
function FBSmets    = getFBScomposition(file)
    if nargin < 1
        FBSmets = {'glucose[s]'
                    'linoleate[s]'
                    'histidine[s]'
                    'isoleucine[s]'
                    'leucine[s]'
                    'lysine[s]'
                    'methionine[s]'
                    'phenylalanine[s]'
                    'threonine[s]'
                    'tryptophan[s]'
                    'valine[s]'
                    'H2O[s]'
                    'O2[s]'
                    'CO2[s]'
                    'alanine[s]'
                    'asparagine[s]'
                    'glutamine[s]'
                    'tyrosine[s]'
                    'cysteine[s]'
                    'arginine[s]'
                    'glycine[s]'
                    'proline[s]'
                    'serine[s]'
                    'aspartate[s]'
                    'glutamate[s]'
                    'Pi[s]'
                    'NH3[s]'
                    'sulfate[s]'
                    'urate[s]'
                    'Fe2+[s]'
                    'Na+[s]'
                    'HCO3-[s]'
                    'H+[s]'
                    'zinc[s]'
                    'K+[s]'
                    'Ca2+[s]'
                    'choline[s]'
                    'ethanolamine[s]'
                    'glycerol[s]'
                    'ornithine[s]'
                    'sulfite[s]'
                    '4-aminobutyrate[s]'
                    'Fe3+[s]'
                    'lactose[s]'
                    'biotin[s]'
                    'acetoacetate[s]'
                    'pyruvate[s]'
                    'L-lactate[s]'
                    'nicotinate[s]'
                    'pantothenate[s]'
                    'folate[s]'
                    'iodide[s]'
                    'nitrite[s]'
                    'chloride[s]'
                    'alpha-tocopherol[s]'
                    'alpha-tocotrienol[s]'
                    'gamma-tocopherol[s]'
                    'gamma-tocotrienol[s]'
                    'hemoglobin[s]'
                    'ascorbate[s]'
                    'thiamin[s]'
                    'homoserine[s]'
                    'malonate[s]'
                    'oxalate[s]'
                    'citrulline[s]'
                    'D-3-amino-isobutanoate[s]'
                    'adenine[s]'
                    'adenosine[s]'
                    'AKG[s]'
                    'AMP[s]'
                    'bilirubin-bisglucuronoside[s]'
                    'bilirubin-monoglucuronoside[s]'
                    'bilirubin[s]'
                    'glycocholate[s]'
                    'taurocholate[s]'
                    'glycochenodeoxycholate[s]'
                    'taurochenodeoxycholate[s]'
                    'cholesterol[s]'
                    'citrate[s]'
                    'CMP[s]'
                    'creatine[s]'
                    'L-carnitine[s]'
                    'cortisol[s]'
                    'cytidine[s]'
                    'deoxycytidine[s]'
                    'deoxyadenosine[s]'
                    'deoxyuridine[s]'
                    'betaine[s]'
                    'GMP[s]'
                    'hypoxanthine[s]'
                    'IMP[s]'
                    'inosine[s]'
                    'nicotinamide[s]'
                    'progesterone[s]'
                    'prostaglandin D2[s]'
                    'prostaglandin E1[s]'
                    'prostaglandin E2[s]'
                    'serotonin[s]'
                    'succinate[s]'
                    'sucrose[s]'
                    'taurine[s]'
                    'thymine[s]'
                    'thymidine[s]'
                    'thyroxine[s]'
                    'triiodothyronine[s]'
                    'testosterone[s]'
                    'UMP[s]'
                    'uracil[s]'
                    'uridine[s]'
                    'creatinine[s]'
                    'D-glucitol[s]'
                    'orotate[s]'
                    'spermine[s]'
                    '2-phospho-D-glycerate[c]'
                    '4-hydroxyphenylpyruvate[c]'
                    '3-hydroxyanthranilate[c]'
                    '3-phospho-D-glycerate[c]'
                    '4-hydroxybenzoate[c]'
                    'trans-4-hydroxy-L-proline[c]'
                    '5-hydroxyindoleacetate[c]'
                    'SAH[c]'
                    'anthranilate[c]'
                    '3-ureidopropionate[c]'
                    'carnosine[c]'
                    'phosphocholine[c]'
                    'L-cystathionine[c]'
                    'dCMP[c]'
                    'DHAP[s]'
                    'dimethylglycine[c]'
                    'ethanolamine-phosphate[c]'
                    'fumarate[c]'
                    'sn-glycerol-3-PC[c]'
                    'glucuronate[c]'
                    'glyceraldehyde[c]'
                    'sn-glycerol-3-phosphate[c]'
                    'glutathione episulfonium ion[c]'
                    'guanidinoacetate[c]'
                    'homocysteine[c]'
                    '3-hydroxy-L-kynurenine[c]'
                    'isocitrate[c]'
                    '4-hydroxy-2-quinolinecarboxylic acid[c]'
                    'L-2-aminoadipate[c]'
                    'kynurenine[c]'
                    'malate[c]'
                    'nicotinate ribonucleotide[c]'
                    'orotidine-5-phosphate[c]'
                    'PEP[c]'
                    'propanoate[s]'
                    'quinolinate[c]'
                    'spermidine[s]'
                    'UDP-glucuronate[c]'
                    'xanthine[c]'
                    'xanthosine-5-phosphate[c]'
                    'xanthosine[c]'
                    'NG,NG-dimethyl-L-arginine[c]'
                    'allantoin[c]'
                    'cis-aconitate[c]'
                    'hippurate[c]'
                    'methylmalonate[c]'
                    'lithocholate[s]'
                    'taurolithocholate[s]'
                    'oleate[s]'
                    'cholesterol-ester pool[c]'};
    else
        fid  = fopen(file,'r');
        scan = textscan(fid,'%s','Delimiter',';');
        fclose(fid);
        FBSmets = scan{1};
    end
end
function outModel   = addCompartmentalizedMetNames(inModel)
% addCompartmentalizedMetNames
% This function adds to the model structure a field called metNamesC
% identical to metNames but with the compartment appended at the end.
%
% INPUT
% inModel:  a genome-scale metabolic model structure. Note: if a field
%           called metNamesC is already in the model, it will be
%           overwritten.
% 
% OUTPUT
% outModel:  a genome-scale metabolic model structure with a metNamesC
%            field.
%
% Usage:  outModel   = addCompartmentalizedMetNames(inModel)
% 2013-07-17 Francesco Gatto
    outModel = inModel;
    if isfield(outModel,'metNamesC')
        outModel = rmfield(outModel,'metNamesC');
    end
    mets                = outModel.metNames;
    nMets               = length(mets);
    comps               = outModel.comps(outModel.metComps);
    metNamesC           = cellfun(@(a,b,c,d) [a,b,c,d],mets,repmat({'['},nMets,1),comps,repmat({']'},nMets,1),'uni',false);
    outModel.metNamesC  = metNamesC;
end
function [metList,metListInd] = getMetListMatchedToUptakeReactions(model,uptakeRxns)    
% function getUptakeReactionsMatchedToMetList
%
% This function scans a model to look for exchange reactions matched to a
% user-defined list of metabolites.
%
% INPUT
% model:    a genome-scale metabolic model. Must be simplified.
% uptakeRxns:   a cell array of rxn IDs or a numeric array of rxn indexes
%               that represent the exchange rxns in the model for which the 
%               matched metabolite is to be searched.
%
% OUTPUT
% metList:  cell array of metabolite names associated to the input list of
%           uptake rxns.
% metListInd:  cell array of metabolite indexes associated to the input list 
%              of uptake rxns.
%
% Usage: metList = getMetListMatchedToUptakeReactions(model,uptakeRxns)
%
% 2013-07-03, Francesco Gatto
    if isfield(model,'unconstrained')
        throw(MException('','Model is still unconstrained, that means metabolites may be in the [x] compartment. Simplify model first.'));
    end
    if ~iscell(uptakeRxns)
        if ~isnumeric(uptakeRxns)
            throw(MException('','Uptake rxns in the input must be either a cell array of rxn IDs or a numeric array of rxn indexes'));
        else
            uptakeRxns = model.rxns(uptakeRxns);
        end
    end
    uptakeRxnsInRefModel = getExchangeRxns(model,'in');
    if ~all(ismember(uptakeRxns,uptakeRxnsInRefModel))
        absentRxnIDsinRefModel = uptakeRxns(~ismember(uptakeRxns,uptakeRxnsInRefModel));
        throw(MException('','Uptake rxn not found in the model: %s\n',absentRxnIDsinRefModel{:}));
    end
    
    nRxns      = numel(uptakeRxns);
    metList    = cell(nRxns,1);
    metListInd = zeros(nRxns,1);
    for i = 1:nRxns
        currentRxn = uptakeRxns(i);
        [~,currentRxnInd] = ismember(currentRxn,model.rxns);    
        metIndex          = find(model.S(:,currentRxnInd));
        if numel(metIndex) < 1            
            throw(MException('','Some uptake rxn IDs in the input have no metabolite associated. Check model'));
        end
        if numel(metIndex) > 1            
            throw(MException('','Some uptake rxn IDs in the input are more than 1 metabolite associated. Check model'));
        else
            metListInd(i) = metIndex;
            metList(i)    = {strcat(model.metNames{metIndex},'[',model.comps{model.metComps(metIndex)},']')};
        end
    end
end
function [metList,metExRxns,modelExRxns] = getUptakeReactionsMatchedToMetList(model,metList,suppressWarnings,suppressOutput)    
% function getUptakeReactionsMatchedToMetList
%
% This function scans a model to look for uptake reactions matched to a
% user-defined list of metabolites. If not found or dead end in the
% extracellular compartment, the function will search whether the
% metabolite is present in the cytosol and add a new uptake reactions.
%
% INPUT
% model:    a genome-scale model. The model must be simplified and contain
%           a field that describes metabolite names together with their
%           compartment ('metNamesC').
% metList:  either a cell array of metabolites or the filename of a plain
%           text file containing metabolite names in each row. Note that
%           metabolites may be either compartmentalized (e.g. glucose[s])
%           or not. In the former scenario, the matched uptake rxn will try
%           to match the corresponding compartment, if such rxn is in the
%           model, or add it to the corresponding compartment, if missing.
% suppressWarnings: a logical. True if no warnings should be displayed.
% suppressOutput:   a logical. True if no output should be displayed.
%
%
% OUTPUT
% metList:  cell array of metabolite names for which an exchange reaction was
%           either found or added to the model
% metExRxns:    cell array of rxn names for the exchange of each metabolite
%               listed in metList
% modelExRxns:  a genome-scale model containing all exchange rxns matched
%               to metList
%
% Usage: [metList,metExRxns,modelExRxns] = getUptakeReactionsMatchedToMetList(model,metList,suppressWarnings)
%
% 2013-06-24, Francesco Gatto
    if nargin < 4 || isempty(suppressOutput)
        suppressOutput = true;
    end
    if nargin < 3 || isempty(suppressWarnings)
        suppressWarnings = false;
    end
    if isfield(model,'unconstrained')
        throw(MException('','Model is still unconstrained, that means metabolites may be in the [x] compartment. Simplify model first.'));
    end
    if ~isfield(model,'metNamesC')
        throw(MException('','There is no field in the model where metabolite names are compartmentalized. Thus metabolites in the list cannot be match. Add metabolite names with compartment first.'));
    end
    if numel(model.metNamesC) ~= numel(model.mets)
        throw(MException('','The field in the model with compartmentalized metabolite names does not match the number of metabolites. Make sure to add this vector after the model has been simplified.'));
    end
    if ~iscell(metList)
        if ischar(metList)
            fid  = fopen(metList,'r');
            scan = textscan(fid,'%s','Delimiter',';');
            metList  = scan{1};
            if isempty(metList)
                throw(MException('','Input text file is expected to be semi-comma delimited.'));
            end
            fclose(fid);
        else
            throw(MException('','Input list of metabolites is expected to be either a cell array of metabolites or a filename for a text file with an array of metabolites.'));
        end
    else
        if all(cellfun(@(a) strcmp(a(end),']'),metList)) && all(cellfun(@(a) strcmp(a(end-2),'['),metList))
            isCompart = true;
        elseif ~any(cellfun(@(a) strcmp(a(end),']'),metList)) || ~any(cellfun(@(a) strcmp(a(end-2),'['),metList))
            isCompart = false;
        else
            throw(MException('','Input list of metabolites is expected to be either all compartmentalized (e.g. glucose[s]) or any of them'));
        end            
    end
    if isCompart
        y = ismember(metList,model.metNamesC);
        if ~all(y)
            throw(MException('','Some metabolites in the input list could not be matched in the current model. Note: check the compartments or remove them'));
        end
    else
        y = ismember(metList,model.metNames);
        if ~all(y)
            throw(MException('','Some metabolites in the input list could not be matched in the current model'));
        end
    end
    uptakeRxns = getExchangeRxns(model,'in'); 
    nMets      = numel(metList);
    metExRxns  = cell(nMets,1);
    toDelete   = [];
    toRename   = [];
    toAdd      = [];
    for i = 1:nMets
        currentMet = metList(i);
        if isCompart
            [~,currentMetInd] = ismember(currentMet,model.metNamesC);
        else
            currentMet = strcat(currentMet,'[s]');
            [~,currentMetInd] = ismember(currentMet,model.metNamesC);
        end
        if currentMetInd ~= 0
            rxnIndexesWmet    = find(model.S(currentMetInd,:));
        else
            rxnIndexesWmet = [];
        end
        if isempty(rxnIndexesWmet)
            % Search for the metabolite in the cytosol
            cytoMet  = {strcat(currentMet{1}(1:end-3),'[c]')};
            [presentInCytosol,currentMetInd] = ismember(cytoMet,model.metNamesC);    
            % If still empty, flag it as to be deleted
            if ~presentInCytosol
                toDelete = [toDelete,i];
                continue
            end
            rxnIndexesWmet    = find(model.S(currentMetInd,:));
            % If still empty, flag it as to be deleted
            if isempty(rxnIndexesWmet)
                toDelete = [toDelete,i];
                continue
            else
                metList(i) = cytoMet;
                toRename = [toRename,i];
            end
        else
            metList(i) = currentMet;
        end
        rxnWMet   = model.rxns(rxnIndexesWmet);
        exRxnWMet = intersect(uptakeRxns,rxnWMet);
        if numel(exRxnWMet) ~= 1
            if numel(exRxnWMet) == 0
                [model,addedRxn]   = addExchangeRxns(model,'in',currentMetInd);
                metExRxns(i) = addedRxn;
                toAdd = [toAdd,i];
            else
                throw(MException('','The metabolite %s has more than one exchange reaction in the model',metList{i}));
            end
        else
            metExRxns(i) = exRxnWMet;
        end
    end
    if ~isempty(toRename) && ~suppressWarnings
        fprintf('\n%s\n','Warning: these metabolites were changed in cytosolic as they were dead end in the extracellular compartment:')
        disp(metList(toRename))
    end
    if ~isempty(toAdd) && ~suppressWarnings
        fprintf('\n%s\n','Warning: these metabolites did not have a matched uptake reactions in the current model, which was therefore added:')
        disp(metList(toAdd))
    end
    if ~isempty(toDelete)
        if ~suppressWarnings
            fprintf('\n%s\n','Warning: these metabolites will be ignored for successive simulations as they do not take part in the network:')
            disp(metList(toDelete))
        end
        metList(toDelete)   = [];
        metExRxns(toDelete) = [];
    end
    if isempty(toRename) && isempty(toAdd) && isempty(toDelete) && ~suppressOutput
        fprintf('\t\t%s\n','All metabolites were succesfully matched to a pre-existing uptake rxn in the model')
    end
    modelExRxns = model;
end
function modelOut = constrainUptRxnsToFluxes(model,metUptRxns,bounds,relaxationTerm)
% function constrainUptRxnsToFluxes
%
% This function constrains the model to open for uptake only those rxns
% listed in metUptRxns. If bounds is provided, an upper limit (and lower
% limit, if bounds is nx2, where n is the number of uptake rxns in the input)
% is imposed to each uptake rxn in the input.
%
% INPUT
% model:        a genome-scale model. Must be simplified.
% metUptRxns:   a nx1 cell array of n uptake rxns contained in model.
% bounds:       a nx1 or nx2 numerical vector array of equality bounds for
%               each uptake rxns in metUptRxns. If nx2 is provided, then
%               the first column is assumed to be a vector of lower bounds,
%               and the second as upper bounds (default is -Inf to +Inf).
%               Note: if n == 1, the same bounds are assumed for all input
%               rxns.
% relaxationTerm: a numeric, indicating the % of relaxation to be used if
%                 "bounds" is a nx1 numeric array (def: 10*eps). Ignored if
%                 "bounds" is nx2.
% OUTPUT
% modelOut:     a genome-scale model in which only those rxns in the input
%               are open for uptake, with the limits provided by bounds (if
%               any).
%
% Francesco Gatto, 2013-06-24
    if nargin < 4
        relaxationTerm = [];
    end     
    if nargin < 3 || isempty(bounds)
        lb = -Inf;
        ub = +Inf;
    else
        if ~isnumeric(bounds)
            throw(MException('','Bounds must be either a nx1 or nx2 numerical vector array'));
        end
        if size(bounds,1)~=numel(metUptRxns)
            if size(bounds,1) ~= 1
                throw(MException('','Bounds must be either a nx1 or nx2 numerical vector array, where n is the number of uptake rxns to constrain'));
            else
                bounds = repmat(bounds,numel(metUptRxns),1);
            end
        end
        if size(bounds,2) == 2
            lb = bounds(:,1);
            ub = bounds(:,2);
        elseif size(bounds,2) == 1
            if isnumeric(relaxationTerm)
                lb = bounds - abs(bounds*relaxationTerm/100);
                ub = bounds + abs(bounds*relaxationTerm/100);
            else
                lb = bounds-10*eps;
                ub = bounds+10*eps;
            end
        else
            throw(MException('','Bounds must be either a nx1 or nx2 vector array'));
        end
    end
    if ~all(ub>=lb)
        throw(MException('','Some upper bounds are not greater or equal to matched lower bounds'));
    end
    uptakeRxns  = getExchangeRxns(model,'in');
    secretRxns  = getExchangeRxns(model,'out');
    model       = setParam(model,'lb',secretRxns,0);
    model       = setParam(model,'ub',uptakeRxns,0);
    model       = setParam(model,'lb',metUptRxns,lb);
    modelOut    = setParam(model,'ub',metUptRxns,ub);
end
function [optim,modelOpt] = optimizeObjective(model,obj,bounds,minmax,minFlux,suppressOutput)
% function optimizeObjective
%
% This function is a wrap-up for solveLP, in that the user can define the
% objective functions and its bounds. Default is maximization.
%
% INPUT
% model:    a genome-scale model
% obj:      a string with the rxn ID of the objective function;
% bounds:   a vector array containing the upper and lower bound for the
%           objective function (default: lb/ub for the rxn as in the model)
% minmax:   a numerical, where 1 is maximization (default) and -1 is
%           minimization.
% minFlux:      determines if a second optimization should be performed
%                 in order to get rid of loops in the flux distribution
%                 0: no such optimization is performed
%                 1: the sum of abs(fluxes) is minimized. This is the
%                 fastest way of getting rid of loops (default)
%                 2: the square of fluxes is minimized. This tends to
%                 distribute fluxes across iso-enzymes, which results in a
%                 larger number of reactions being used
%                 3: the number of fluxes is minimized. This can result
%                 in the flux distributions that are the easiest to
%                 interpret. Note that this optimization can be very slow
%                 (opt, default 1)
% suppressOutput: a logical indicated if output on screen should be
%           suppressed (def: true).
%
% OUTPUT
% optim:    a structure containing the solution of the optimization:
%               x: the solution vector of fluxes;
%               f: the value of the objective functuion;
%               stat: a flag for the optimization problem;
% modelOpt: the model with settings for the objective function implemented.
%
% Usage: [optim,modelOpt] = optimizeObjective(model,obj,bounds,minmax)
%
% 2013-10-22 Francesco Gatto
    if ~ismember(obj,model.rxns)
        throw(MException('','Objective function (rxnID) not found in the model.'));
    end
    if nargin < 6 || isempty(suppressOutput)
        suppressOutput = true;
    end
    if nargin < 5 || isempty(minFlux)
        minFlux = 1;
    end
    if nargin < 4 || isempty(minmax)
        minmax = 1;
    end
    if nargin < 3 || isempty(bounds)
        objInd = getIndexes(model,obj,'rxns');
        bounds = [model.lb(objInd),model.ub(objInd)];
    end
    modelOpt = setParam(model,'obj',obj,minmax);
    modelOpt = setParam(modelOpt,'lb',obj,bounds(1));
    modelOpt = setParam(modelOpt,'ub',obj,bounds(2));
    optim    = solveLP(modelOpt,minFlux);    
    if ~suppressOutput    
        if optim.stat == 1
            fprintf('\t%s\n','Solution for the given objective with the given constraints: Found!')
            fprintf('\t\t%s\t%d%s%6.4f%s%6.4f%s\n','The objective function optimum was:',optim.x(getIndexes(modelOpt,obj,'rxns')),' [',bounds(1),'-',bounds(2),' ]')
            fprintf('\t\t%s\n','Using the following exchange fluxes:')
            printFluxesTidy(modelOpt,optim.x,true,[],[],'\t\t\t%rxnID (%rxnName):%flux [%lower %upper]\n')
        else
            fprintf('\t%s\n','Solution for the given objective with the given constraints: Not found!') 
        end
    end 
end
function printFluxesTidy(model, fluxes, onlyExchange, cutOffFlux, outputFile,outputString,metaboliteList)
if nargin<3
    onlyExchange=true;
end
if nargin<4
    cutOffFlux=10^-8;
end
if isempty(cutOffFlux)
    cutOffFlux=10^-8;    
end
if nargin<5
    fid=1;
else
    if ~isempty(outputFile)
        fid=fopen(outputFile,'w');
    else
        fid=1;
    end
end
if nargin<6
    outputString='%rxnID (%rxnName):%flux\n';
end
if isempty(outputString)
    outputString='%rxnID (%rxnName):%flux\n';
end
if nargin<7
    metaboliteList={};
end
if numel(fluxes)~=numel(model.rxns)
   throw(MException('','The number of fluxes and the number of reactions must be the same.')); 
end

%Only keep reactions involving the defined metabolites
if ~isempty(metaboliteList)
    I=ismember(upper(model.metNames),upper(metaboliteList));
    [crap K]=find(model.S(I,:));
    
    %Delete all other reactions
    toDelete=true(numel(model.rxns),1);
    toDelete(K)=false;
    model=removeRxns(model,toDelete);
    fluxes(toDelete)=[];
end

%Remove reactions which are below the cut off
toDelete=abs(fluxes)<cutOffFlux;
model=removeReactions(model,toDelete,true,true);
fluxes(toDelete)=[];

if any(strfind(outputString,'%eqn'))
    %Construct the equations
    eqn=constructEquations(model);
else
    eqn=cell(numel(model.rxns),1);
    eqn(:)={''};
end
if any(strfind(outputString,'%element'))
    %For printing equations using the composition
    cModel=model;
    cModel.metNames=cModel.metFormulas;
    cModel.metNames(cellfun(@isempty,cModel.metNames))={'?'};
    element=constructEquations(cModel);
else
    element=cell(numel(model.rxns),1);
    element(:)={''};
end

if any(strfind(outputString,'%unbalanced')) || any(strfind(outputString,'%lumped'))
   balanceStructure=getElementalBalance(model);
end

unbalanced=cell(numel(model.rxns),1);
unbalanced(:)={''};
if any(strfind(outputString,'%unbalanced'))
   unbalanced(balanceStructure.balanceStatus==0)={'(*)'};
   unbalanced(balanceStructure.balanceStatus<0)={'(-)'};
end

lumped=cell(numel(model.rxns),1);
lumped(:)={''};
if any(strfind(outputString,'%lumped'))
    compounds={'C' 'P' 'S' 'N' 'R' 'O' 'H'};
    for i=1:numel(model.rxns)
        leftGroup='';
        rightGroup='';
        for j=1:7
            I=balanceStructure.leftComp(i,j);
            if I~=0
                if I==1
                    leftGroup=[leftGroup compounds{j}];
                else
                    leftGroup=[leftGroup compounds{j} num2str(I)];
                end
            end
            I=balanceStructure.rightComp(i,j);
            if I~=0
                if I==1
                    rightGroup=[rightGroup compounds{j}];
                else
                    rightGroup=[rightGroup compounds{j} num2str(I)];
                end
            end
        end
        if model.rev(i)
        	lumped{i}=[leftGroup ' <=> ' rightGroup];
        else
            lumped{i}=[leftGroup ' => ' rightGroup]; 
        end
    end
end

for i=1:numel(model.rxns)
   %Only print if it's an exchange reaction or if all reactions should be
   %printed. Exchange reactions only have reactants or only products.
   reactants=model.S(:,i)<0;
   products=model.S(:,i)>0;
   
   %Only print if the absolute value is >= cutOffFlux
   if (onlyExchange==false || (~any(reactants) || ~any(products)))
       printString=outputString;
        
       %Produce the final string
       printString=strrep(printString,'%rxnID',model.rxns{i});
       printString=strrep(printString,'%eqn',eqn{i});
       printString=strrep(printString,'%rxnName',model.rxnNames{i});
       printString=strrep(printString,'%lower',num2str(model.lb(i)));
       printString=strrep(printString,'%upper',num2str(model.ub(i)));
       printString=strrep(printString,'%obj',num2str(model.c(i)));
       printString=strrep(printString,'%flux',num2str(fluxes(i)));
       printString=strrep(printString,'%element',element{i});
       printString=strrep(printString,'%unbalanced',unbalanced{i});
       printString=strrep(printString,'%lumped',lumped{i});
       fprintf(fid,printString);
   end
end

if fid~=1
    fprintf('File successfully saved.\n');
    fclose(fid);
end
end
function geneEssentiality = applyGenePerturbation(model,perturbation,settings)
% applyGenePerturbation
%
% This function checks the effect of a single or double gene perturbation
% on the objective function of a model. It outputs a structure containing 
% whether the perturbation was ineffective, redundant, toxic, or lethal. It
% features some advanced options to exclude some genes from the analysis.
% This may also be used to reduce the list of double gene perturbation to
% test.
%
% INPUT
% model:    a genome-scale metabolic model structure.
% perturbation: a string indicating which perturbation to apply. Either
%               'sgd' or 'dgd'.
% settings: a structure containing the following fields:
%           - genePerturbationMethod: a string indicating whether 'fba' or
%                                 'moma' should be used when applying 'sgd'
%                                 or 'dgd' (def: 'fba').
%           - singleGenesToExclude: a nx1 cell array of strings indicating which genes in the
%                             model should be excluded from being perturbed (def:
%                             none).
%           - doubleGenesToExclude: a nx2 cell array of strings indicating which gene pairs in the
%                             model should be excluded from being perturbed (def:
%                             none). Only if "perturbation" is 'dgd'.
%           - excludeSingleGenesNotToPerturbFromDgd: a logical indicating
%                             whether the list of genes to pair to test DGD
%                             should be obtained by all the genes in the
%                             model except those indicated in
%                             "singleGenesNotToPerturb" (def: false). Only 
%                             if "perturbation" is 'dgd'.
%           - suppressOutput: a logical if output on screen should be
%                             suppressed (def: true).
%           - saveSGDfilename: a string if SGD output should be saved on
%                              file with this name (def: '').
%
% OUTPUT
% geneEssentiality: a structure containing the following fields:
%   - redundantGenes: a nx1 cell array of strings
%       containing all genes that cannot affect the matched rxns since iso-
%       enzymes are reported in the model for the same rxn.
%   - deadEndGenes: a nx1 cell array of strings
%       containing all genes that cannot affect the matched rxns since the
%       rxn(s) cannot carry flux (i.e. dead-ends).
%   - excludedGenes: a nx1 cell array of strings
%       containing all genes that were not perturbed (identical to input
%       argument "singleGenesToExclude").
%   - viableGenes/viableDoubleGenes: a nx1 (/nx2) cell array of strings
%       containing all genes (/double combinations of genes) that if deleted
%       do not affect the optimization of the model objective function.
%   - lethalGenes/lethalSingleGenes: a nx1 cell array of strings
%       containing all genes that if deleted abolish the optimization of 
%       the model objective function. 
%   - lethalDoubleGenes: a nx2 cell array of strings containing all double 
%       combinations of genes that if deleted abolish the optimization of 
%       the model objective function. Only if 'dgd' was selected.
%   - toxicGenes/toxicDoubleGenes: a nx1 (/nx2) cell array of strings 
%       containing all genes (/double combinations of genes) that if deleted 
%       diminish the optimal value of the model objective function.
%
% NOTE#1: This function disregards complexes. Any one gene can encode a
%           reaction even if parts of the complex is deleted. 
% NOTE#2: If "excludeSingleGenesNotToPerturbFromDgd" is false, first not
%           lethal single gene deletions are computed and then double gene
%           deletions starting from this list are tested. If true, the
%           first step is skipped, and the list of gene pairs to test for
%           double gene deletions is obtained directly from all genes
%           except those contained in "singleGenesNotToPerturb". If
%           possible, supply a meaningful list of "singleGenesNotToPerturb"
%           to save computational time.
%
% Usage: geneEssentiality = applyGenePerturbation(model,perturbation,settings)
%
% 2014-08-26 Francesco Gatto
    if nargin < 3
        genePerturbationMethod  = 'fba';
        excludeSingleGenesNotToPerturbFromDgd = false;
        doubleGenesToExclude    = [];
        singleGenesToExclude    = [];
        suppressOutput          = true;
    else
        if isstruct(settings)
            if isfield(settings,'genePerturbationMethod')
                genePerturbationMethod = settings.genePerturbationMethod;
                if ~strcmp(genePerturbationMethod,'fba') && ~strcmp(genePerturbationMethod,'moma')
                    throw(MException('','Gene perturbation method must either be "fba" or "moma"'));
                end
            else
                genePerturbationMethod  = 'fba';
            end                
            if isfield(settings,'excludeSingleGenesNotToPerturbFromDgd')
                excludeSingleGenesNotToPerturbFromDgd = settings.excludeSingleGenesNotToPerturbFromDgd;
                if ~islogical(excludeSingleGenesNotToPerturbFromDgd)
                    throw(MException('','excludeSingleGenesNotToPerturbFromDgd must be a cell array of strings.'));
                end
            else
                excludeSingleGenesNotToPerturbFromDgd = false;
            end            
            if isfield(settings,'doubleGenesToExclude')
                doubleGenesToExclude = settings.doubleGenesToExclude;
                if ~isempty(doubleGenesToExclude)
                    if ~iscell(doubleGenesToExclude)                          
                        throw(MException('','Genes to exclude must be a cell array of strings.'));
                    end
                    if size(doubleGenesToExclude,2) ~= 2
                        throw(MException('','Double genes to exclude must be a nx2 cell array of strings.'));
                    else
                        firstGene  = doubleGenesToExclude(:,1);
                        secondGene = doubleGenesToExclude(:,2);
                        if ~all(ismember(firstGene,model.genes)) || ~all(ismember(secondGene,model.genes))
                            throw(MException('','Some of the double genes to exclude were not found in the model.'));
                        end
                    end       
                else
                    firstGene  = [];
                    secondGene = []; 
                end
            else
                doubleGenesToExclude = [];
                firstGene  = [];
                secondGene = []; 
            end
            if isfield(settings,'singleGenesToExclude')
                singleGenesToExclude = settings.singleGenesToExclude;
                if ~isempty(singleGenesToExclude)
                    if ~iscell(singleGenesToExclude)                          
                        throw(MException('','Genes to exclude must be a cell array of strings.'));
                    end
                    if ~all(ismember(singleGenesToExclude,model.genes))
                        throw(MException('','Some genes to exclude were not found in the model.'));
                    end
                end
            else
                singleGenesToExclude = [];
            end 
            if isfield(settings,'suppressOutput')
                suppressOutput = settings.suppressOutput;
                if ~isempty(suppressOutput)
                    if ~islogical(suppressOutput)                          
                        throw(MException('','suppressOutput must be a logical.'));
                    end
                end
            else
                suppressOutput = true;
            end 
            if isfield(settings,'saveSGDfilename')
                saveSGDfilename = settings.saveSGDfilename;
                if ~isempty(saveSGDfilename)
                    if ~ischar(saveSGDfilename)                          
                        throw(MException('','saveSGDfilename must be a string.'));
                    end
                end
            else
                saveSGDfilename = '';
            end 
        else
            disp('Warning: "settings" is not a structure. All contents are ignored.')
            genePerturbationMethod  = 'fba';
            excludeSingleGenesNotToPerturbFromDgd = false;
            doubleGenesToExclude    = [];
            singleGenesToExclude    = [];
            suppressOutput          = true;
            saveSGDfilename         = '';
        end
    end
    if nargin < 2
        throw(MException('','"perturbation" must be provided as an input.'));
    else
        if ~ischar(perturbation)
            throw(MException('','"perturbation" must be a string'));
        else
            if ~any(strcmp(perturbation,{'sgd';'dgd'}))
                throw(MException('','"perturbation" must be either "sgd" or "dgd"'));
            end
        end
    end
    if isempty(model.c)
        throw(MException('','The model does not have any objective function implemented.'));
    end
%   testGeneDeletions works on indexes rather than gene list. Convert gene
%   names to gene indexes
    singleGenesToExcludeInd = getIndexes(model,singleGenesToExclude,'genes');
    doubleGenesToExcludeInd = [getIndexes(model,firstGene,'genes'),getIndexes(model,secondGene,'genes')];
    [genes fluxes originalGenes details]= testGeneDeletions(model,perturbation,genePerturbationMethod,singleGenesToExcludeInd,doubleGenesToExcludeInd,...
        excludeSingleGenesNotToPerturbFromDgd,suppressOutput,saveSGDfilename);
    if strcmp(perturbation,'sgd')
        notEssentialGenes = genes;
        toxicGenes        = getToxicGenes(model,originalGenes,notEssentialGenes,fluxes);
        geneEssentiality.viableGenes       = setdiff(originalGenes(notEssentialGenes),toxicGenes);
        geneEssentiality.lethalGenes       = originalGenes(details==2);
        geneEssentiality.toxicGenes        = toxicGenes;
        geneEssentiality.redundantGenes    = originalGenes(details==3);
        geneEssentiality.deadEndGenes      = originalGenes(details==4);
        geneEssentiality.excludedGenes     = originalGenes(singleGenesToExcludeInd);
    else
        geneEssentiality.redundantGenes      = originalGenes(details==3);
        geneEssentiality.deadEndGenes        = originalGenes(details==4);
        geneEssentiality.excludedGenes       = setdiff(originalGenes(singleGenesToExcludeInd),originalGenes(details==3|details==4));
        geneEssentiality.excludedDoubleGenes = doubleGenesToExclude;
        geneEssentiality.lethalSingleGenes   = originalGenes(details==2);
        geneEssentiality.lethalDoubleGenes   = originalGenes(genes(sum(fluxes,1)==0,:));
        geneEssentiality.viableDoubleGenes   = originalGenes(genes(fluxes(model.c~=0,:)==model.ub(model.c~=0),:));
        geneEssentiality.toxicDoubleGenes    = originalGenes(genes(sum(fluxes,1)~=0&fluxes(model.c~=0,:)<0.5*model.ub(model.c~=0),:));
    end
end
function [genes fluxes originalGenes details]= testGeneDeletions(model,testType,analysisType,essentialSingleGeneInds,essentialDoubleGeneInds,excludeSingleGenesNotToPerturbFromDgd,suppressOutput,saveSGDfilename)
% testGeneDeletions
%   It wraps up the function findGeneDeletions in RAVEN by adding the
%   possibility to exclude some genes from the analysis
    originalModel=model;
    if strcmpi(analysisType,'moma')
        refModel = model;
    end
    originalGenes=model.genes;
    details=zeros(numel(model.genes),1);
    %First simplify the model to reduce the size
    model=simplifyModel(model,true,false,true,true);
    model=removeReactions(model,{},true,true); %Removes unused genes
    details(~ismember(originalGenes,model.genes))=4;
    [~, geneMapping]=ismember(model.genes,originalGenes);
    %Get the genes that should be deleted. Not all genes are deleted in all
    %optimizations. For example, if all reactions involving a gene have
    %iso-enzymes a single deletion of that gene will never have an effect.
    %FBA and SGD:   All genes that encode reactions on their own
    %FBA and DGD:   All genes that have at most one other gene encoding for at
    %               least one of its reactions and doesn't participate
    %               exclusively in reactions where a single deletion results in
    %               an unsolvable problem
    %MOMA and SGD:  All genes that encode reactions on their own
    %MOMA and DGD:  All combinations of genes that have at most one other gene 
    %               encoding for at east one of its reactions and doesn't 
    %               participate exclusively in reactions where a single 
    %               deletion results in an unsolvable problem
    genesForSGD=[];
    if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
        if strcmpi(testType,'sgd')
            I=sum(model.rxnGeneMat,2)>1; %Reactions with iso-enzymes
        else
            I=sum(model.rxnGeneMat,2)>2; %Reactions with more than one iso-enzyme
        end

        %Find genes involved only in these reactions
        [~, inIsoRxns]=find(model.rxnGeneMat(I,:));
        [~, inNonIsoRxns]=find(model.rxnGeneMat(~I,:));
        genesForSGD=unique(inNonIsoRxns);
        details(geneMapping(setdiff(inIsoRxns,inNonIsoRxns)))=3;
    end
    %Get the genes that should be deleted in a SGD
    if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
       if ~isempty(essentialSingleGeneInds)
           [I,J] = ismember(essentialSingleGeneInds,geneMapping(genesForSGD));
           genesToModify=genesForSGD(~ismember(geneMapping(genesForSGD),essentialSingleGeneInds));
           details(essentialSingleGeneInds(I))=5;
       else
           genesToModify=genesForSGD; 
       end
    end
    %Output statistics so far
    if ~suppressOutput
        fprintf('\t%s\t%s\n\t%s\t%s\n','Perturbation type:',testType,'Perturbation method:',analysisType)
        fprintf('\t%s\t%d\n','N° genes in the model:',numel(originalGenes))
        fprintf('\t%s\t%d\n','N° dead-end genes:',sum(details==4))
        fprintf('\t%s\t%d\n','N° redundant genes:',sum(details==3))
        fprintf('\t%s\t%d\n','N° excluded genes:',sum(details==5))
        fprintf('\t%s\t%d\n','Total number of genes to be perturbed:',numel(genesToModify))
    end
    %Start file if output should be saved
    if strcmpi(testType,'sgd') && ~isempty(saveSGDfilename)
        fid = fopen(saveSGDfilename,'wt');
    end
    %Do single deletion/over expression. This is done here since the double
    %deletion depends on which single deletions prove lethal
    %(to reduce the size of the system). Note that if a list of genes to
    %pair for dgd is provided, this step is skipped.
    if strcmpi(testType,'sgd') || (strcmpi(testType,'dgd') && ~excludeSingleGenesNotToPerturbFromDgd)
        if ~suppressOutput
            fprintf('\t\t%s\n','Running single gene perturbation test...')
        end
        fluxes=zeros(numel(model.rxns),numel(genesToModify));
        solvable=true(numel(genesToModify),1);
        for i=1:numel(genesToModify)
           if ~suppressOutput
                fprintf('\t\t\t%s%d%s\t%s%s','Testing single gene perturbation #',i,':',model.genes{genesToModify(i)},'... ')
           end
           if ~isempty(saveSGDfilename)
               fprintf(fid,'%s\t',model.genes{genesToModify(i)});
           end
           if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
               %Constrain all reactions univocally encoded by the gene to 0
               tempModel=applyGeneKO(model,genesToModify(i),true);
           else
               %To over express a gene, the stoichiometry of the corresponding
               %reactions are changed so that the same flux leads to a higher
               %production
               tempModel=model;
               tempModel.S(:,I)=tempModel.S(:,I).*oeFactor;
           end
           if strcmpi(analysisType,'fba') || strcmpi(testType,'dgd')
                sol=solveLP(tempModel);
           else
                [fluxA , ~, flag]=qMOMA(tempModel,refModel);
                sol.x=fluxA;
                sol.stat=flag;
           end
           %If the optimization terminated successfully
           if sol.stat==1
               fluxes(:,i)=sol.x;
               details(geneMapping(genesToModify(i)))=1;
               if ~suppressOutput
                    fprintf('%s\n','Viable or toxic!')
               end
               if ~isempty(saveSGDfilename)
                   fprintf(fid,'%4.2f\n',-sol.f);
               end
           else
               solvable(i)=false;
               details(geneMapping(genesToModify(i)))=2;
               if ~suppressOutput
                    fprintf('%s\n','Lethal!')
               end
               if ~isempty(saveSGDfilename)
                   fprintf(fid,'%f\n',NaN);
               end
           end
        end
        fluxes=fluxes(:,solvable);
        genes=geneMapping(genesToModify(solvable));
        if ~suppressOutput
            fprintf('\t\t%s\n','Single gene perturbation complete!')
            fprintf('\t\t%s\t%d\n','N° viable (or toxic) genes:',sum(details==1))
            fprintf('\t\t%s\t%d\n','N° lethal genes:',sum(details==2))
        end
    end
    %For double deletions FBA or MOMA
    if strcmpi(testType,'dgd')
        %This is a little lazy but it's fine. Check which genes that have
        %already been labled as either ony with too many iso-enzymes or
        %non-solveable as a single deletion.
        if ~excludeSingleGenesNotToPerturbFromDgd
            genesForDGD = originalGenes(details==1);
            if ~suppressOutput
                fprintf('\t\t%s\n','Running double gene perturbation test on previously viable genes...')
                fprintf('\t\t%s\t%d\n','Total number of genes to be perturbed:',numel(genesForDGD))
            end
        else
            % In this case, the list of double gene deletions to test is
            % derived directly from the modifiable list of genes available
            % for SGD, which should be properly filtered using
            % "singleGenesToExclude"
            genesForDGD = originalGenes(geneMapping(genesToModify));
            if ~suppressOutput
                fprintf('\t\t%s\n','Running double gene perturbation test...')
            end
        end
        [~, I] = ismember(genesForDGD,model.genes);
        genesToModify = nchoosek(I,2);
        if ~suppressOutput
            fprintf('\t\t%s\t%d\n','Total number of gene pairs:',size(genesToModify,1))
        end
        if ~isempty(essentialDoubleGeneInds)
            genesToModify(ismember(geneMapping(genesToModify),essentialDoubleGeneInds,'rows'),:) = [];
            if ~suppressOutput
                fprintf('\t\t%s\t%d\n','Total number of gene pairs excluded by users:',numel(ismember(geneMapping(genesToModify),essentialDoubleGeneInds,'rows')))
                fprintf('\t\t%s\t%d\n','Total number of gene pairs to be perturbed:',numel(genesToModify))
            end
        end
        genes=geneMapping(genesToModify);
        fluxes=sparse(numel(model.rxns),size(genesToModify,1));
        for i=1:size(genesToModify,1)
           [I crap]=find(model.rxnGeneMat(:,genesToModify(i,:)));
           if ~suppressOutput
                fprintf('\t\t\t%s%d%s\t%s%s','Testing double gene perturbation #',i,':',model.genes{genesToModify(i,:)},'... ')
           end
          %Constrain all reactions involving the gene to 0
           tempModel=setParam(model,'eq',model.rxns(I),0);

           if strcmpi(analysisType,'fba')
                sol=solveLP(tempModel);
           else
                [fluxA , ~, flag]=qMOMA(tempModel,refModel);
                sol.x=fluxA;
                sol.stat=flag;
           end

           if sol.stat==1
               fluxes(:,i)=sol.x;
               if ~suppressOutput
                    fprintf('%s\n','Viable or toxic!')
               end
           elseif ~suppressOutput
               fprintf('%s\n','Lethal!')
           end
        end
        if ~suppressOutput
            fprintf('\t\t%s\n','Double gene perturbation complete!')    
        end
    end
    %Map back to the old model
    [~, I]=ismember(model.rxns,originalModel.rxns);
    temp=fluxes;
    fluxes=sparse(numel(originalModel.rxns),size(temp,2));
    fluxes(I,:)=temp;
end
function modelKO = applyGeneKO(model,genes,onlyUnivocally)
% applyGeneKO
%
% Returns a model containing as many gene-knockouts as in genes. If
% onlyUnivocally is set to false, then all rxns encoded by the genes are
% bounded to zero (else, only those that are univocally encoded by
% the genes)
%
% model: a genome-scale model structure:
% genes: a cell array of strings containing the genes to knockout.
% onlyUnivocally: a logical. If true, only rxns that are univocally encoded 
%                 by the genes are bounded to 0 (opt, def: true).
% modelKO: a genome-scale model structure with all genes in genes knocked
%          out
% Usage: modelKO = applyGeneKO(model,genes)
%
% Francesco Gatto, 2013-11-22
if nargin < 3
    onlyUnivocally = true;
end
   nGenes = numel(genes);
   for i = 1:nGenes
       currentGene = genes(i);
       try
           geneInds = getIndexes(model,currentGene,'genes');
       catch  %#ok<CTCH>
           throw(MException('','One or more genes is not in the model.'));
       end
       I       = find(model.rxnGeneMat(:,geneInds));
       I       = unique(I);
       if onlyUnivocally
           u       = sum(model.rxnGeneMat(I,:),2)==1;
           model   = setParam(model,'eq',model.rxns(I(u)),0);
       else
           model   = setParam(model,'eq',model.rxns(I),0);
       end
   end
   modelKO = model;
end
function toxicGenes = getToxicGenes(model,originalGenes,notEssentialGenes,fluxes,cutoff)
    objInd = find(model.c);
    objVal = fluxes(objInd,:);
    objUb  = model.ub(objInd);
    objLb  = model.lb(objInd);
    if nargin < 3
        objUb = cutoff;
    end
    % Definition of toxic genes is those genes that are not
    % lethal (i.e. the objective has not even reached the lower bound) but
    % are at most half the maximum upper bound. Numerical errors produce
    % effectively 0 flux, but sometimes lower than 0. Therefore it is safer
    % only to include the limit on the upper bound and check for
    % truthfulness later
    toxicGenes = originalGenes(notEssentialGenes(objVal < (objUb / 2)));
end
function [toxicDataset,geneEssentiality] = verifyToxic(modelOpt,obj,objBounds,geneEssentiality,toxicThreshold)
 % Verify nature of toxic genes
     if ~isempty(geneEssentiality.toxicGenes)
        solWT = optimizeObjective(modelOpt,obj,objBounds,1,false);
        maxGrowthWT = solWT.f;
        for i = 1:numel(geneEssentiality.toxicGenes)
            toxicG  = geneEssentiality.toxicGenes(i);
            modelKO = applyGeneKO(modelOpt,toxicG);
            sol     = optimizeObjective(modelKO,obj,objBounds,1,false);
            maxGrowth(i) = sol.f/maxGrowthWT; %#ok<*NASGU>
        end
        toxicDataset = dataset(geneEssentiality.toxicGenes,maxGrowth','VarNames',{'Genes','MaxGrowth'});
        geneEssentiality.toxicGenes  = toxicDataset.Genes(toxicDataset.MaxGrowth<toxicThreshold);
        geneEssentiality.viableGenes = union(geneEssentiality.viableGenes,toxicDataset.Genes(toxicDataset.MaxGrowth>=toxicThreshold));
     else
        toxicDataset = dataset([]);
     end
end
function h = drawPerturbationSummary(essentialityStructure,perturbation,figurename)
% drawPerturbationSummary
%
% This function draw a pie chart summarizing the overall effect of each
% perturbation (either 'sgd', 'dgd', or 'am') in a genome-scale metabolic
% model.
%
% INPUT
% essentialityStructure: a perturbation essentiality structure as output
%                        from applyGenePerturbation or
%                        applyMetabolitePerturbation.
% perturbation: a string containing either 'sgd', 'dgd', or 'am' (def:
%               empty, the function guesses the perturbation type)
% figurename: a string or cell array of 
%
% OUTPUT
% h: a handle for the pie chart
%
% Usage: 
% 
% 2013-07-11 Francesco Gatto
    if nargin < 3
       figurename = [];
    else
        if ~isempty(figurename) && ~ischar(figurename)    
            if ~iscell(figurename)
                throw(MException('','If a figure name is provided, a valid string must be provided'));
            else
                figurename = figurename{:};
                if ~ischar(figurename)
                    throw(MException('','If a figure name is provided, a valid cell array of strings must be provided'));
                end
            end 
        end
    end
    if nargin < 2 || isempty(perturbation)
        if isfield(essentialityStructure,'lethalGenes')
            perturbation = 'sgd';
        elseif isfield(essentialityStructure,'lethalSingleGenes')
            perturbation = 'dgd';
        elseif isfield(essentialityStructure,'lethalMetabolite')
            perturbation = 'am';
        else
            throw(MException('','Essentiality structure appears not to contain all fields to draw a summary'));
        end
    end
    switch perturbation
        case 'sgd'
            nViable  = numel(essentialityStructure.viableGenes);
            nToxic   = numel(essentialityStructure.toxicGenes);
            nLethal  = numel(essentialityStructure.lethalGenes);
            nRedund  = numel(essentialityStructure.redundantGenes);
            nDeadEnd = numel(essentialityStructure.deadEndGenes);
            nExclude = numel(essentialityStructure.excludedGenes);
            explode  = [0,0,0,0,1,1];
            pie([nRedund,nDeadEnd,nExclude,nViable,nToxic,nLethal],explode,{'Redundant';'DeadEnds';'Excluded';'Viable';'Toxic';'Lethal'})
        case 'dgd'
            nViable         = numel(essentialityStructure.viableDoubleGenes);
            nToxic          = numel(essentialityStructure.toxicDoubleGenes);
            nLethalSingle   = numel(essentialityStructure.lethalSingleGenes);
            nLethalDouble   = numel(essentialityStructure.lethalDoubleGenes);
            nRedund         = numel(essentialityStructure.redundantGenes);
            nDeadEnd        = numel(essentialityStructure.deadEndGenes);
            nExclude        = numel(essentialityStructure.excludedGenes);
            explode  = [0,0,0,1,1,1];
            pie([nRedund,nDeadEnd,nExclude,nViable,nToxic,nLethalSingle,nLethalDouble],explode,{'Redundant';'DeadEnds';'Excluded';'Viable';'Toxic';'LethalSingle';'LethalDouble'})
        case 'am'
            return
    end
    if ~isempty(figurename)
        print('-painters','-dpdf','-r300',figurename)
    end
end