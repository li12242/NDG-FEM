function results = runStdCellTest()
%> @brief Run all the unit tests of the NdgCell module.
%>
%> The suite covers the StdCell classes (StdCellTest), the pure-Matlab
%> Polylib gateway ports (PolylibGatewayTest) and the triangular prism cell
%> (StdPrismTriTest). Paths are resolved relative to this file, so the
%> suite runs from any current folder.
%>
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

thisDir = fileparts(mfilename('fullpath'));
srcRoot = fullfile(thisDir, '..', '..'); % source/NDG-FEM-master
addpath( fullfile(srcRoot, 'lib') );
addpath( fullfile(srcRoot, 'NdgCell') );

import matlab.unittest.TestSuite
suite = [ TestSuite.fromFile( fullfile(thisDir, 'StdCellTest.m') ), ...
          TestSuite.fromFile( fullfile(thisDir, 'PolylibGatewayTest.m') ), ...
          TestSuite.fromFile( fullfile(thisDir, 'StdPrismTriTest.m') ) ];
results = suite.run();
end% func
