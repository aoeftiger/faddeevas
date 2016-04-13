for x = 0:0.5:20
for y = 0:0.5:20
  z=fadf(x+y*i);
  printf("fadf(%.2f, %.2f) = %.17e, %.17e\n", x, y, real(z), imag(z));
endfor
printf("\n")
endfor

%y=0:0.5:10;
%for x = 0:0.5:10
%  z=fadf(x+y*i);
%  for i = 1:length(y)
%    printf("fadf(%.2f, %.2f) = %.17e, %.17e\n", x, y(i), real(z(i)), imag(z(i)));
%  endfor
% printf("\n")
%endfor