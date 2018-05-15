function val = issameset(A,B)
% function bool  = issameset(A,B);
% compares if vectors A and B contain the same set of numbers
% regardless of ordering

val = 0;

% if lengths are not equal, then it wont work
if length(A) ~= length(B)
    val = 0;
    return;
end

% test basic in order 
if A == B
    val = 1;
    return
end

% ok, now loop over A and compare with elements of B
for i=1:length(A)
   if find (A(i) == B)
       % so far so good.
       val = 1;
       continue;
   else
       % if the element isnt there, then return false
       val = 0;
       return;
   end
end


