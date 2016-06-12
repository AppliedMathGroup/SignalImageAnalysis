% MATLAB function for the analysis of time series, using the Square Wave Method
% Developed by: Ricardo E. Monge, Sherry Gapper, Osvaldo Skliar
% Date: June 6, 2016 
% The MATLAB command described in [1] performs a time series analysis using 
% the Square Wave Transform, and is accessible from the MATLAB interface 
% under the name swt.
%
% swt is a user-defined function requiring three parameters to function correctly. 
% The first parameter V corresponds to the time series to be analyzed, input as a 
% standard MATLAB vector. The second parameter f corresponds to the sampling
% frequency in Hz of the data in parameter V. Finally, the third parameter Dt 
% corresponds to the time interval in seconds of the entire data set. 
%
% The outputs are a two-column matrix S, and a single numerical value Dm,
% indicating the approximation quality (i.e., the maximum of the moduli of the
% differences between the original data values and the data values computed
% by adding the respective square waves). S will have a row for each data value
% in the input data set, in which the value in the first column corresponds to
% the square wave frequency (Fi in [1]), and the value in the second column
% corresponds to the square wave coefficient (Ci in [1]).
%
% References
% [1] Skliar O., Monge R. E., Oviedo G. and Gapper S. (2016). A New Method for the 
%     Analysis of Signals: The Square Wave Transform, Revista de Matemática. 
%     Teoría y aplicaciones 2016, Vol. 23(1), pp. 85-110.
%
function [S,Dm]=swt(V,f,Dt) 
 n=length(V);
 C=zeros(n);
 %  CREATE THE COMPUTATION MATRIX
 for i=1:n, 
  for j=1:n,
   lj=n-j+1;
   q=idivide(int32(i),int32(lj));
   r=mod(i,lj);
   qmod2=mod(q,2);
   if (qmod2==0 && r==0),
    C(i,j)=-1;
   elseif (qmod2==0 && r~=0), 
    C(i,j)=1;
   elseif (qmod2==1 && r==0),
    C(i,j)=1;
   elseif (qmod2==1 && r~=0), 
    C(i,j)=-1;
   end 
  end
 end
 % SOLVE THE SYSTEM OF EQUATIONS
 SV=linsolve(C,V');
 for i=1:n,
  SF(i,j)=1/(2*Dt) * (n/(n-i+1));
 end
 % COMPUTE APPROXIMATION
 SEQUENCE=[1:n];
 SQWAVELENGTHS=2*Dt*SEQUENCE/n;
 CV=zeros(1,n);
 de=Dt/n;
 for i=1:n,
  ti=i*de+de/2;
  for j=1:n,
    CV(i)=CV(i)+ SV(i)*power(-1,floor(2*ti/SQWAVELENGTHS(j)));
  end  
 end
 % COMPUTE DM
 Dm=max(abs(V-CV));
 % PREPARE S
 S=[SF; SV];
end
