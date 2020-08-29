function [outputArg1,outputArg2] = QuadraticProgram(inputArg1,inputArg2)
% Solves min l*cx + 1/2x'Cx, s.t. Ax=b, x>=0 and l is a scalar
% exitflag is 0 if solved successfully, 1 if infeasible and -1 if unbounded
% Performs a Phase 1 procedure...
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

function [] = ShortFormSolve()

%w = b;
%cdual = -l*p;
%id = -l*p>=0;
%zp = cdual(id);
%zn = cduals(~id);
Ai = [A eye(m,m) zeros(m,2*n)];
Ai = [Ai; ([C zeros(n, m) ones(n, n) -ones(n, n)])];
bi = [b; -l*p];
ci = [zeros(n,1); ones(m,1); zeros(n,1); zeros(n,1)];
[~, xi, pi, indices, exitflag] = RSM(Ai, bi, ci, m, n);

x = xi(indices([1:n,m+n:m+3*n]));


% Initial xb and basic cost
cb=ci(indices([1:n,m+n:m+3*n]));
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

