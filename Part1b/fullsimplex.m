function [z, x, pi, indices, exitflag] = fullsimplex(A, b, c, m, n)
% Solves min cx s.t. Ax=b, x>=0
% exitflag is 0 if solved successfully, 1 if infeasible and -1 if unbounded
% Performs a Phase I procedure starting with an all artificial basis
% using simplexPh, and then uses it for phase 2

% Phase 1
Ph=1;
% Initialise
indices=(n+[1:m])'; % All artificial basis
Bmatrix=eye(m);
ct=zeros(n+m,1);
ct(n+1:n+m)=1; % Artificial costs
[z, x, pi, indices, ~, Bmatrix] = simplexPH(A, b, ct, m, n, Bmatrix, indices, Ph);

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
    [z, x, pi, indices, exitflag, ~] = simplexPH(A, b, ct2, m, n, Bmatrix, indices, Ph);
    
    % If degenrate artificials, don't show artificials in final solution
    if (Ph==2)&&(any(indices>n))&&(exitflag==0)
        x=x(1:n);
    end
end
end

function [z, x, pi, indices, exitflag, Bmatrix] = simplexPH(A, b, c, ~, n, Bmatrix, indices, Ph)
% Solves min cx s.t. Ax=b, x>=0, using the Phase 1,2 approach
% starting with basic variables listed in vector indices
% and basis matrix Bmatrix
% exitflag is 0 if solved successfully and -1 if unbounded
% returns optimal vector x and its value z, along with pi, 
% indices of basic variables, and the final basis matrix

% Initialise
cb=c(indices);

% Simplex step until an exit condition is met
while true
    
    % Initialise and update basic solution
    IBmatrix=inv(Bmatrix);
    xb=IBmatrix*b;
    
    % Entering step
    pi=(cb'*IBmatrix)';
    % Don't consider artificials
    [as, cs, s]=findenter(A,pi,c(1:n));
    if s==0 % If at optimal sol, don't preform leaving step, end
        exitflag=0;
        x=zeros(n,1); % Initialise and assign the optimal solution
        x(indices)=xb;
        z=cb'*xb;
        
        return
    end
    
    % Leaving step
    leave = findleaveEXT(Bmatrix, as, xb, n, indices, Ph);
    if leave==0 % If the problem is unbounded, end simplex
        exitflag=-1;
        
        % Avoid non assignment issues
        x=zeros(n,1); % Initialise and assign the end solution
        x(indices)=xb;
        z=cb'*xb;
        return
    end
    
    % Initialise for the next simplex step
    [Bmatrix, indices, cb] = update(Bmatrix, indices, cb, cs, as, s, leave);
end
end

function [as, cs, s] = findenter(Amatrix, pi, c)
% Given the complete m by n matrix Amatrix,
% the complete cost vector c with n components
% the vector pi with m components
% findenter finds the index of the entering variable and its column
% It returns the column as, its cost coefficient cs, and index s
% Returns s=0 if no entering variable can be found (i.e. optimal)

% Find the minimum reduced cost
[m,s]=min(c'-pi'*Amatrix);
% If this minimum is not negative, the solution is optimal
% NOTE: I have included a tolerance, to deal with near-singular Bmatrices
% and Matlabs precision. This fixes the example included in BFSStudentv1
if m>=-1.0e-6
    s=0;
    as=0; % Avoid non assignment issues
    cs=0;

else % Grab the column and cost of the entering var
    as=Amatrix(:,s);
    cs=c(s);
end
end

function [leave] = findleaveEXT(Bmatrix, as, xb, n, indices, Ph)
% Given entering column as and vector xb of basic variables
% findleaveEXT finds a leaving column of basis matrix Bmatrix
% It preferentially selects artificial vars if in Phase 2
% It returns leave=0 if no column can be found (i.e. unbounded)

% Compute the denominator vector of the ratio test
bottom=(Bmatrix\as)';

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

function [newBmatrix, newindices, newcb] = update(Bmatrix, indices, cb, cs, as, s, leave)
% Bmatrix is current m by m basis matrix
% indices is a column vector current identifiers for basic variables in
% order of B columns
% cb is a column vector of basic costs in the order of B columns
% as is the entering column
% s is the index of the entering variable
% leave is the column index (p) of the basis matrix that must leave
% (not its variable index t)
% update replaces column leave of Bmatrix with as to give newBmatrix
% replaces row leave of indices with enter to give newindices
% replaces row leave of cb with cs to give newcb

newBmatrix=Bmatrix;
newBmatrix(:,leave)=as;

newindices=indices;
newindices(leave)=s;

newcb=cb;
newcb(leave)=cs;
end