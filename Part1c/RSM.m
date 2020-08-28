function [z, x, pi, indices, exitflag] = RSM(A, b, c, m, n)
% Solves min cx s.t. Ax=b, x>=0
% exitflag is 0 if solved successfully, 1 if infeasible and -1 if unbounded
% Performs a Phase 1 procedure starting with an all artificial basis
% by calling rsmSTP, and then uses it for Phase 2

% Phase 1
Ph=1;
% Initialise
indices=(n+[1:m])'; % All artificial basis
IBmatrix=eye(m);
ct=zeros(n+m,1);
ct(n+1:n+m)=1; % Artificial costs
[z, x, pi, indices, ~, IBmatrix] = rsmSTEP(A, b, ct, m, n, IBmatrix, indices, Ph);

% Check for problems before starting Phase 2
if z>0
    exitflag=1;
    return
  
else % Phase 2
    Ph=2;
    if any(indices>n)
        ct2=c;
        ct2(indices(indices>n))=0; % Artificial costs
    else
        ct2=c;
    end
    [z, x, pi, indices, exitflag, ~] = rsmSTEP(A, b, ct2, m, n, IBmatrix, indices, Ph);
    
    % If degenerate artificials, don't show artificials in final solution
    if (Ph==2)&&(any(indices>n))&&(exitflag==0)
        x=x(1:n);
    end
end
end

function [z, x, pi, indices, exitflag, IBmatrix] = rsmSTEP(A, b, c, ~, n, IBmatrix, indices, Ph)
% Solves min cx s.t. Ax=b, x>=0
% starting with basic variables listed in vector indices
% and basis matrix Bmatrix
% exitflag is 0 if solved successfully and -1 if unbounded
% returns optimal vector x and its value z, along with pi, 
% indices of basic variables and the final inverse basis matrix

% Initial xb and basic cost
cb=c(indices);
xb=IBmatrix*b;

% Revised simplex step until an exit condition is met
while true
    
    % Entering step
    pi=(cb'*IBmatrix)';
    % Don't consider artificials
    [as, cs, s]=findenterIP(A,pi,c(1:n),indices(indices<=n));
    if s==0 % If at optimal sol, don't preform leaving step, end simplex
        exitflag=0;
        x=zeros(n,1); % Initialise and assign the optimal solution
        x(indices)=xb;
        z=cb'*xb;
        
        return
    end
    
    % Leaving step
    leave = rsmfindleaveEXT(IBmatrix, as, xb, n, indices, Ph);
    if leave==0 % If the problem is unbounded, end simplex
        exitflag=-1;

        % Avoid non assignment issues
        x=zeros(n,1); % Initialise and assign the end solution
        x(indices)=xb;
        z=cb'*xb;
        return
    end
    
    % Initialise for the next simplex step
    [IBmatrix, indices, cb, xb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave, xb);
end
end

function [as, cs, s] = findenterIP(Amatrix, pi, c, indices)
% Given the complete m by n matrix Amatrix,
% the complete cost vector c with n components
% the vector pi with m components
% findenterIP finds the index of the entering variable and its column,
% using individual pricing, ignoring basic columns
% It returns the column as, its cost coefficient cs, and index s
% Returns s=0 if no entering variable can be found (i.e. optimal)

if length(c)>length(indices) % Check there is a non basic col to price
    AcolI=1:length(c);
    for i=AcolI(~ismember(AcolI,indices)) % Only price non basics
        
        % Find the reduced cost
        rc=c(i)-pi'*Amatrix(:,i);
        
        % If this cost is negative, choose it to enter
        % NOTE: tolerance for machine zero
        if rc<-1.0e-6
            
            % Grab the column and cost of the entering var
            as=Amatrix(:,i);
            cs=c(i);
            s=i;
            return
        end
    end
end

% If no non basic vars have a negative reduced cost, the sol is optimal 
s=0;
as=0; % Avoid non assignment issues
cs=0;
end

function [leave] = rsmfindleaveEXT(IBmatrix, as, xb, n, indices, Ph)
% Given entering column as and vector xb of basic variables
% rsmfindleaveEXT finds a leaving column of the original basis matrix, Bmatrix
% It preferentially selects artificial vars if in Phase 2
% It returns leave=0 if no column can be found (i.e. unbounded)

% Compute the denominator vector of the ratio test
bottom=(IBmatrix*as)';

% Extended algorithm
if Ph==2&&(any(indices>n))
    indexes=1:n;
    for i=indexes(indices>n) % Cycle through artificial indexes
        if (bottom(i)~=0)&&(indices(i)>n)
            leave=i;
            return
        end
    end
end    

% Normal findleave 

% Make any non valid denominators not get chosen for the minimum index
xb(bottom<=0)=nan;

% (B^-1)b is just xb in std form. Acquire the minimum ratio's index
[~,mindex]=min((xb)'./bottom);

% Check if vars can be decreased without bound
if all(bottom<=0)
    leave=0;
else
    leave=mindex;
end
end

function [newIBmatrix, newindices, newcb, newxb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave, xb)
% IBmatrix is current m by m inverse basis matrix
% indices is a column vector current identifiers for basic variables in
% order of B columns
% cb is a column vector of basic costs in the order of B columns
% as is the entering column
% s is the index of the entering variable
% xb is the old basic solution
% leave is the column index (p) of the basis matrix that must leave
% (not its variable index t)
% updateGJ uses gauss-jordan pivoting to update the inverse basis with as,
% giving newIBmatrix and newxb
% replaces row leave of indices with enter to give newindices
% replaces row leave of cb with cs to give newcb

% Preform row reduction on the system to get the new inverse
Augm=[xb IBmatrix IBmatrix*as];

redRow=Augm(leave,:)/Augm(leave,end); % Find the multiples of the reducer
Augm=Augm-Augm(:,end)*redRow; % Row reduce the system so final column is 0
Augm(leave,:)=redRow; % Add the row used for this back in
newIBmatrix=Augm(:,2:end-1);

newxb=Augm(:,1);

newindices=indices;
newindices(leave)=s;

newcb=cb;
newcb(leave)=cs;

end

%{
function [as, cs, s] = findenterPP(Amatrix, pi, c, indices)
% Given the complete m by n matrix Amatrix,
% the complete cost vector c with n components
% the vector pi with m components
% findenterPP finds the index of the entering variable and its column,
% using partial pricing, breaking the matrix into N blocks,
% and ignoring basic columns
% It returns the column as, its cost coefficient cs, and index s
% Returns s=0 if no entering variable can be found (i.e. optimal)

Amatrix(:,indices)=[]; % Ignore basic cols

if isempty(Amatrix) % All A columns are basic, optimal
    s=0;
    as=0; % Avoid non assignment issues
    cs=0;
    return
end

c(indices)=[];

n=size(Amatrix,2);

% PLEASE NOTE: I originally wrote this for partial pricing becuase this
% wasn't specified in the documentation, but since we must do individual
% pricing I used findenterPP
N=25; % The number of blocks to price
blCols=ceil(n/N); % Number columns in each block

% Cycle through all but the last block/semi-block
for i=0:(max(floor(n/blCols)-2,0))
    
    % Find the minimum reduced cost in this block
    [mi,s]=min(c(i*blCols+1:(i+1)*blCols)'-pi'*Amatrix(:,i*blCols+1:(i+1)*blCols));
    
    % If this minimum is negative, choose it as the entering var
    if mi<-1.0e-6 % NOTE: tolerance to deal with near-singular matrices
        s=s+i*blCols; % Index of the minimum wrt all blocks
        as=Amatrix(:,s); % Grab the column and cost of the entering var
        cs=c(s);
        
        Ci=1:(n+length(indices)); % Represent all original columns that shift
        Ci(indices)=[]; % Remove shifted columns 
        s=Ci(s); % Index wrt to all original columns
        return
    end
    % If the minimum is not negative, try the next block
end

% Check the last block/semi-block to find its minimum reduced cost
[mi,s]=min(c((i+1)*blCols+1:end)'-pi'*Amatrix(:,(i+1)*blCols+1:end));

% If this minimum is negative, choose it as the entering var
if mi<-1.0e-6 % NOTE: tolerance to deal with near-singular matrices
    s=s+(i+1)*blCols; % Index of the minimum wrt all blocks
    as=Amatrix(:,s); % Grab the column and cost of the entering var
    cs=c(s);
    
    Ci=1:(n+length(indices)); % Represent all original columns that shift
    Ci(indices)=[]; % Remove shifted columns 
    s=Ci(s); % Index wrt to all original columns
    
else % Since there are no negative reduced costs, the solution is optimal
    s=0;
    as=0; % Avoid non assignment issues
    cs=0;
end
end
%}