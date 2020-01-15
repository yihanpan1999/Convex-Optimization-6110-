function X = myInpaint(Xh, Omega)

%%% Figure out relevant size (n, m) from the two input arguments.
n = length(Xh);
m = length(Omega);

%%% Construct the E matrix and S matrix.
temp = speye(n-1);
left = [temp sparse(n-1,1)];
right = [sparse(n-1,1) -temp];
D = left + right;
up = kron(speye(n),D);
down = kron(D,speye(n));
E = [up; down];
full = speye(n^2);
S = full(Omega,:);

%%% Construct the A matrix and f and b vectors.
identity = speye(2*n*(n-1));
zero = sparse(m,2*n*(n-1));
A = [E, -identity, identity;
    S, zero, zero];
f = [sparse(1,n^2), ones(1,4*n^2-4*n)];
x = Xh(:);
b = [sparse(2*n*(n-1),1);x(Omega,:)];

%%% Call MATLAB linprog to solve the linear program to get z.
lb = sparse(5*n^2-4*n,1);
optimset('LargeScale','on');
options = optimoptions('linprog','Algorithm','interior-point','Display','none','ConstraintTolerance',1e-2,'OptimalityTolerance',1e-2,'MaxIterations',10);
[z,object] = linprog(f,[],[],A,b,lb,[],options);

%%% Extract the image x from z and shape it into a matrix X for output.
X = reshape(z(1:n^2,:),n,n);

