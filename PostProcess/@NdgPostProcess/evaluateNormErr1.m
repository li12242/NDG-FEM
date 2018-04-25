%======================================================================
%> @brief Calculate the L2-norm error.
%>
%> The L1-norm error is calculated from the integral formula
%>
%> \f$ L^2(c) = \frac{1}{A}\sqrt{ \sum_{k=1}^K \int_{\Omega_k}
%> { \left( u_h - u_{ext} \right) } } \f$
%> (Bunya, Kubatko and Westerink, et al., COMPUT METHOD APPL M, 2009)
%> 
%> where exact solution should be stored in the properties `fext`, and  
%> the error are calculated from the formula for each variable field.
%>
%> @retval err Vector of the L2 norm error for each variable field.
%======================================================================
%> This function is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ err ] = evaluateNormErr1( obj, fphys, fext )
err = zeros(obj.Nvar, 1);

for m = 1:obj.Nmesh
    totalArea = sum( obj.meshUnion(m).LAV );
    for fld = 1:obj.Nvar
        temp = fphys{m}(:,:,fld) - fext{m}(:,:,fld);
        absErr = sum( abs( temp(:) ) );
        err(fld) = err(fld) + absErr/totalArea;
    end
end
end% func