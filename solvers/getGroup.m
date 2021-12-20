function [newA,G,ind,reorder] = getGroup(A,grpNUM)
% A,b      - data
% grpNUM   - number of groups

% G,ind    - input for Def_P.m(weight wi = sqrt(|G_i|))
[p] = size(A,2); %p:number of variables
if grpNUM > p, error('grpNUM > #variables'); end
ind = zeros(3,grpNUM);
G = [1:p];

%% reordering p variables
reorder = randperm(p);
newA = A(:,reorder);

%% partition [1:p] into grpNUM groups

% probs = rand(nG,1); probs = probs / sum(probs); 
% probs = ones(nG,1)/nG;

probs = 1 + 0.3 * sign(randn(grpNUM,1)) .* rand(grpNUM,1); 
probs = probs / sum(probs); 
probs = cumsum(probs);

for i = 1:grpNUM
    if i == 1
        tmp = round(probs(1)*p); tmp = max(tmp,1);
        ind(1,1) = 1; ind(2,1) = tmp; ind(3,1) = sqrt(tmp);
    else
        ind(1,i) = ind(2,i-1) + 1;
        ind(2,i) = max(round(probs(i)*p),ind(1,i));
        ind(3,i) = sqrt(ind(2,i)-ind(1,i));
    end
end


end
