% analyze_bus_errors.m
%
% Takes solutions on two different W1-line defect domains and reports their
% difference along only the bus itself, disregarding the rods. Assume the
% line defects are 15 blocks wide.
%
% For convenience, define A = SOLUTION1 and B = SOLUTION2 beforehand.

% Get the number of buffer rows for each.
An = ((length(A) / 15) - 1) / 2; Bn = ((length(B) / 15) - 1) / 2;

ERROR = struct('U', cell(1,15));

x = An*15;
y = Bn*15;

% err = zeros(1,15);
for i=1:15
    ERROR(i).U = A(x+i).U - B(y+i).U;
%     err(i) = max(max(ERROR(i).U));
end

% disp(max(err));