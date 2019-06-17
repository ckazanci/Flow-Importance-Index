function [Npossible,NpossibleF,NpossibleFcond,conditions] = ...
    searchFlows(Npossible,NpossibleF,NpossibleFcond,conditions, ...
                S,reducedS,construct,currentRow,currentFlows,numFlows)

numRows = rows(construct);
if((numRows==columns(construct))&&(numRows>0))
    % This is a square matrix. Test to see if this flow is acceptable.
    if(rank(construct)==numRows)
        % This set of flows uniquely determine all of the other
        % flows. We can use this as a case to add to the list of
        % acceptable flow.
        Npossible = Npossible + 1;
        conditions( Npossible ) = cond(construct);
        
        % Now for each flow in this list increment the number of
        % sets of flows that this flow appears in.
        working_flow_set = currentFlows; %setdiff(1:numFlows,currentFlows);
        NpossibleF(working_flow_set) = NpossibleF(working_flow_set) + 1;
        NpossibleFcond(working_flow_set) = NpossibleFcond(working_flow_set) + 1/conditions( Npossible );
    end
else
    % This is not a square matrix, and we need to add another flow
    % to the set before testing whether or not it is an acceptable
    % set of flows.
    [numRows,numCols] = size(reducedS);
    for lupe=1:numCols
        % Go through each column in the RREF. Add the column if it
        % has a non-zero entry.
        if(reducedS(currentRow,lupe)!=0)
            % This column might work. Add it to the list and test
            % it out.
            %printf('Checking column %d for row %d\n',lupe,currentRow)
            currentFlows(currentRow) = lupe;
            [Npossible,NpossibleF,NpossibleFcond,conditions] = ...
                searchFlows( Npossible,NpossibleF,NpossibleFcond,conditions, ...
                             [S(:,1:(lupe-1)) S(:,(lupe+1:numCols))], ...
                             [reducedS(:,1:(lupe-1)) reducedS(:,(lupe+1:numCols))], ...
                             [construct S(:,lupe)], ...
                             currentRow+1, ...
                             currentFlows,numFlows);
        end
    end
end
