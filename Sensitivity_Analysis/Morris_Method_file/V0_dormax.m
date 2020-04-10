function y = V0_dormax(p)

M = 18 ; % number of model parameters
p1 = p(:);
if ~isnumeric(p1); error('input argument ''p_'' must be numeric'); end
if length(p1)~=M; error('input argument ''p_'' must have %d components',M); end

[mineralA,transcriptsA,DNA_totA,DNA_a,DNA_i,DNA_s,RNA,P,CO2,rate,ty,SSE] = sensitivity_V0(p1);

% OUTPUTS

y(1)    = max(mineralA);                    % Maximum mineralization first application
y(2)    = max(transcriptsA);                % Maximum transcripts first application
y(3)    = max(DNA_totA);                    % Maximum total DNA first application
y(4)    = max(RNA);                         % Maximum RNA per gene first application
y(5)    = min(P);                           % Minimum pestidice first application
y(6)    = max(CO2);                         % Maximum CO2 first application
y(7)    = max(DNA_a);                       % Maximum CO2 first application
y(8)    = max(DNA_i);                       % Maximum CO2 first application
y(9)    = max(DNA_s);                       % Maximum CO2 first application
y(10)   = SSE;                              % Maximum CO2 first application
a       = ty([find(RNA==max(RNA))]);        % Time of maximum mineralization rate or degradation rate
if isempty(a)==1
    y(11)   = NaN;
else
    y(11)   = a(1);                         % Time of maximum gene expression and maximum gene abundance
end
b       = ty([find(rate==max(rate))]);     
if isempty(b)==1
    y(12)   = NaN;
else
    y(12)   = b(1);                         % Time of maximum mineralization rate or degradation rate
end
c       = (max(mineralA)/2);     
[val,idx]       = min(abs(mineralA-c));
d       = ty(idx);     
y(13)   = d(1);                             % Time at which 50% of maximum mineralization is reached