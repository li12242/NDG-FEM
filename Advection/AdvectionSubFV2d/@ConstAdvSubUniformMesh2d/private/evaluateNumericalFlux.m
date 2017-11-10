function [ flux ] = evaluateNumericalFlux( fm, fp, nxm, nym, u, v )

[ uNorm ] = u .* nxm + v .* nym;
flux = ( fm.*( sign(uNorm) + 1 )*0.5 + fp.*( 1 - sign(uNorm)  )*0.5 ) .* uNorm;

end

