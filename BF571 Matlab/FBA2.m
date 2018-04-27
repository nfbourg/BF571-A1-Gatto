function [v_sol optimal_val status] = FBA(S , Lbound, Ubound , OBJ , Metabolites , Reactions)
    
% Modification of Daniel Segre's script for FBA using glpkmex
% Further Modification, 2017, to use Matlab linprog
% Developed by Evan Snitkin, Sep. 6, 2006
%
% Takes as input :
%	S 		- The stoichiometric matrix, rows are metabolites and columns are reactions
%	Ubound		- A vector the same length as the number of reactions, indicating the upper bound for each
%	Lbound		- A vector the same length as the number of reactions, indicating the lower bound for each
%	OBJ		- A vector the same length as the number of reactions, indicating which fluxes should be optimized
%	Metabolites	- A vector of metabolite names
%	Reactions	- A vector of reaction names
%   OutputMode - boolean variable ; if 1, print out details of reavctions
%
%	
% Returns : The fluxes through each reaction and the exit status of glpk.
% Do help glpkmex for more details on glpkmex parameters
%
	
    OutputMode=1;
	
	%SET THE B VECTOR FOR THE RIGHT HAND SIDE OF EQUATION, Sv= 0
	Bvector = zeros(size(S,1) , 1);
	
	%NUMBER OF METABOLITES (ROWS OF S, CONSTRAINTS)
	M = size(S,1); 
	
	%NUMBER OF REACTIONS (COLUMNS OF S, VARIABLES)
	N = size(S,2); 	
	
	% -1 MAXIMIZES, 1 MINIMIZES
	min_max = -1;

    % RUN OPTIMIZATION (Matlab linprog)
    options = optimoptions('linprog','Algorithm','dual-simplex');
	[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(min_max*OBJ,[],[],S,Bvector,Lbound,Ubound,options);
    
    v_sol = X;
    optimal_val = min_max*FVAL;
    status = EXITFLAG;
    
    
    

    % ------------- PRETTY OUTPUT OF ALL FLUXES INFO -----------------------------------

    % Generates cell arry Reactions of strings containing details of each reaction
    % substrates and products

    React_Names=Reactions;
    LB=Lbound;
    UB=Ubound;

    for f=1:N % loop on all reactions
       curr_column=S(:,f);
       curr_index=find(curr_column);
       tmpstr=[''];
       for g=1:length(curr_index)
           tmpstr=[tmpstr '[' num2str(curr_column(curr_index(g)))  ' ' Metabolites{curr_index(g)} ']'];
           PrettyReactions{f}=tmpstr;
       end   
    end


%     disp(' ');
%     disp('REACTIONS and BOUNDS')
%     for f=1:length(React_Names)
%         [string2print,error_sprintf]=sprintf('%d\tName=%s\tReaction={%s}\tBounds=[%.2f,%.2f]\tFlux=%.2f',f,React_Names{f},Reactions{f},LB(f),UB(f),v_sol(f));
%         disp(string2print);
%     end



    %EXTRA OUTPUT: For each metabolite, writes all fluxes consuming and producing it:
    disp_metabolites=0; % Display info only if this variable is set to 1
    if (disp_metabolites)
        for f=1:M;
            disp(' ')
            disp(['Metabolite ' num2str(f) ' = ' Metabolites{f}])
            pcf= find(S(f,:));
            curr_f=v_sol(pcf).*abs(S(f,pcf))';
            for g=1:length(curr_f)
                disp([ '    '  num2str(pcf(g)) '  '  React_Names{pcf(g)} '  '  num2str(curr_f(g))])
            end
        end
    end






