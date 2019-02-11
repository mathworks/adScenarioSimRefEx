function currSmartCount = markedSmartActorCount()
% Keep track of number of actors marked as smart.

% Copyright 2018 The MathWorks, Inc.

persistent count;
if isempty(count)
    count = 0;
end
count = count + 1;
currSmartCount = count;
end