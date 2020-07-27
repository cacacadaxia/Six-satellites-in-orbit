%  this function solves Kepler's equation,
%  computing E as a function of M and e
%
function E = EofMe(M,e,tol)
   if ( nargin<3 ), tol=1e-11; end
 En  = M;
 En1 = En - (En-e*sin(En)-M)/(1-e*cos(En));
 while ( abs(En1-En) > tol )
  En = En1;
  En1 = En - (En-e*sin(En)-M)/(1-e*cos(En));
 end;
 E = En1;