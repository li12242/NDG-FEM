function [ physPost ] = makeNdgPostProcessFromNdgPhys( phys )
physPost = NdgPostProcess( phys.meshUnion, phys.getOption('outputNetcdfCaseName') );
end

