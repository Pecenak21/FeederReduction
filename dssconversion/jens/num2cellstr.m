   function C = num2strcell(A)
   % NUM2STRCELL Convert an array to a cell array of strings
   % C = num2strcell(A)
   % e.g., num2strcell(1:3) = { '1', '2', '3' }

   n = length(A);
   C = cell(size(A));
   for i = 1:prod(size(A))
      C{i} = num2str(A(i));
   end