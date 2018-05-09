function [A] = A_factor(nrel)
% This subroutine wants to calculate the factor A that accounts Fresnel
% reflections
% nrel relative refractive index
      if nrel > 1.0
       A=504.332889-2641.00214*nrel+5923.699064*nrel*nrel-7376.355814*nrel*nrel*nrel+5507.5304*nrel^4-2463.357945*nrel^5+610.956547*nrel^6-64.8047*nrel^7;     
      elseif nrel==1.0
       A=1.0; 
      elseif nrel < 1.0
       A=3.084635-6.531194*nrel+8.357854*nrel*nrel-5.082751*nrel*nrel*nrel+1.171382*nrel^4;
      end
end