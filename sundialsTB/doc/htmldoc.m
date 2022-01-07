% HTMLDOC - Creates html documentation for sundialsTB
%

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2022, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2006/07/17 16:49:49 $

% Set location of sundialsTB template files location
s = fileparts(which(mfilename));
tmpl = fullfile(s,'html_files/stb_template');

% If the output directory does not exist, create it
system('mkdir -p stb_guide');

% Set output dir
cd('..');
doc_dir = 'doc/stb_guide';

%-----------------------------
% INSTALL and SETUP scripts
%-----------------------------

top_files = {'./install_STB.m'};

m2html('mfiles',top_files,...
       'recursive', 'off',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','foo');

%-----------------------------
% CVODES
%-----------------------------

cmd = sprintf('rm -f %s/cvodes.html',doc_dir);
system(cmd);

cvodes_fct = {'cvodes/CVodeSetOptions.m'...
              'cvodes/CVodeQuadSetOptions.m'...
              'cvodes/CVodeSensSetOptions.m'...
              'cvodes/CVodeInit.m'...
              'cvodes/CVodeQuadInit.m'...
              'cvodes/CVodeSensInit.m'...
              'cvodes/CVodeAdjInit.m'...
              'cvodes/CVodeInitB.m'...
              'cvodes/CVodeQuadInitB.m'...
              'cvodes/CVodeReInit.m'...
              'cvodes/CVodeQuadReInit.m'...
              'cvodes/CVodeSensReInit.m'...
              'cvodes/CVodeAdjReInit.m'...
              'cvodes/CVodeReInitB.m'...
              'cvodes/CVodeQuadReInitB.m'...
              'cvodes/CVode.m'...
              'cvodes/CVodeB.m'...
              'cvodes/CVodeSensToggleOff.m'...
              'cvodes/CVodeGetStats.m'...
              'cvodes/CVodeGetStatsB.m'...
              'cvodes/CVodeGet.m'...
              'cvodes/CVodeSet.m'...
              'cvodes/CVodeSetB.m'...
              'cvodes/CVodeFree.m'...
              'cvodes/function_types/CVRhsFn.m'...
              'cvodes/function_types/CVSensRhsFn.m'...
              'cvodes/function_types/CVQuadRhsFn.m'...
              'cvodes/function_types/CVRootFn.m'...
              'cvodes/function_types/CVDenseJacFn.m'...
              'cvodes/function_types/CVBandJacFn.m'...
              'cvodes/function_types/CVJacTimesVecFn.m'...
              'cvodes/function_types/CVPrecSetupFn.m'...
              'cvodes/function_types/CVPrecSolveFn.m'...
              'cvodes/function_types/CVGcommFn.m'...
              'cvodes/function_types/CVGlocalFn.m'...
              'cvodes/function_types/CVMonitorFn.m'...
              'cvodes/function_types/CVRhsFnB.m'...
              'cvodes/function_types/CVQuadRhsFnB.m'...
              'cvodes/function_types/CVDenseJacFnB.m'...
              'cvodes/function_types/CVBandJacFnB.m'...
              'cvodes/function_types/CVJacTimesVecFnB.m'...
              'cvodes/function_types/CVPrecSetupFnB.m'...
              'cvodes/function_types/CVPrecSolveFnB.m'...
              'cvodes/function_types/CVGcommFnB.m'...
              'cvodes/function_types/CVGlocalFnB.m'...
              'cvodes/function_types/CVMonitorFnB.m'};

% First run over all of them to get the cross-references right
m2html('mfiles','cvodes',...
       'recursive','on',...
       'globalHypertextLinks','on',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','cvodes');

% Re-run without examples to remove source
m2html('mfiles',cvodes_fct,...
       'htmldir',doc_dir,...
       'source','off', ...
       'template',tmpl,...
       'indexFile','cvodes');


%-----------------------------
% IDAS
%-----------------------------

cmd = sprintf('rm -f %s/idas.html',doc_dir);
system(cmd);

idas_fct = {'idas/IDASetOptions.m'...
            'idas/IDAQuadSetOptions.m'...
            'idas/IDASensSetOptions.m'...
            'idas/IDAInit.m'...
            'idas/IDAQuadInit.m'...
            'idas/IDASensInit.m'...
            'idas/IDAAdjInit.m'...
            'idas/IDAInitB.m'...
            'idas/IDAQuadInitB.m'...
            'idas/IDAReInit.m'...
            'idas/IDAQuadReInit.m'...
            'idas/IDASensReInit.m'...
            'idas/IDAAdjReInit.m'...
            'idas/IDAReInitB.m'...
            'idas/IDAQuadReInitB.m'...
            'idas/IDACalcIC.m'...
            'idas/IDACalcICB.m'...
            'idas/IDASolve.m'...
            'idas/IDASolveB.m'...
            'idas/IDASensToggleOff.m'...
            'idas/IDAGetStats.m'...
            'idas/IDAGetStatsB.m'...
            'idas/IDAGet.m'...
            'idas/IDASet.m'...
            'idas/IDASetB.m'...
            'idas/IDAFree.m'... 
            'idas/function_types/IDAResFn.m'...
            'idas/function_types/IDASensResFn.m'...
            'idas/function_types/IDAQuadRhsFn.m'...
            'idas/function_types/IDARootFn.m'...
            'idas/function_types/IDADenseJacFn.m'...
            'idas/function_types/IDABandJacFn.m'...
            'idas/function_types/IDAJacTimesVecFn.m'...
            'idas/function_types/IDAPrecSetupFn.m'...
            'idas/function_types/IDAPrecSolveFn.m'...
            'idas/function_types/IDAGcommFn.m'...
            'idas/function_types/IDAGlocalFn.m'...
            'idas/function_types/IDAMonitorFn.m'...
            'idas/function_types/IDAResFnB.m'...
            'idas/function_types/IDAQuadRhsFnB.m'...
            'idas/function_types/IDADenseJacFnB.m'...
            'idas/function_types/IDABandJacFnB.m'...
            'idas/function_types/IDAJacTimesVecFnB.m'...
            'idas/function_types/IDAPrecSetupFnB.m'...
            'idas/function_types/IDAPrecSolveFnB.m'...
            'idas/function_types/IDAGcommFnB.m'...
            'idas/function_types/IDAGlocalFnB.m'...
            'idas/function_types/IDAMonitorFnB.m'};

% First run over all of them to get the cross-references right
m2html('mfiles','idas',...
       'recursive','on',...
       'globalHypertextLinks','on',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','idas');

% Re-run without examples to remove source
m2html('mfiles',idas_fct,...
       'htmldir',doc_dir,...
       'source','off', ...
       'template',tmpl,...
       'indexFile','idas');


%-----------------------------
% KINSOL
%-----------------------------

cmd = sprintf('rm -f %s/kinsol.html',doc_dir);
system(cmd);

kinsol_fct = {'kinsol/KINSetOptions.m'...
              'kinsol/KINInit.m'...
              'kinsol/KINSol.m'...
              'kinsol/KINGetStats.m'...
              'kinsol/KINFree.m'...
              'kinsol/function_types/KINSysFn.m'...
              'kinsol/function_types/KINDenseJacFn.m'...
              'kinsol/function_types/KINBandJacFn.m'...
              'kinsol/function_types/KINJacTimesVecFn.m'...
              'kinsol/function_types/KINPrecSetupFn.m'...
              'kinsol/function_types/KINPrecSolveFn.m'...
              'kinsol/function_types/KINGcommFn.m'...
              'kinsol/function_types/KINGlocalFn.m'};
 
% First run over all of them to get the cross-references right
m2html('mfiles','kinsol',...
       'recursive','on',...
       'globalHypertextLinks','on',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','kinsol');

% Re-run without examples to remove source
m2html('mfiles',kinsol_fct,...
       'htmldir',doc_dir,...
       'source','off', ...
       'template',tmpl,...
       'indexFile','kinsol');

%-----------------------------
% NVECTOR
%-----------------------------

cmd = sprintf('rm -f %s/nvector.html',doc_dir);
system(cmd);

nvector_fct = {'nvector/N_VDotProd.m'...
               'nvector/N_VL1Norm.m'...
               'nvector/N_VMax.m'...
               'nvector/N_VMaxNorm.m'...
               'nvector/N_VMin.m'...
               'nvector/N_VWL2Norm.m'...
               'nvector/N_VWrmsNorm.m'};

m2html('mfiles',nvector_fct,...
       'recursive','on',...
       'globalHypertextLinks','on',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','nvector');

%-----------------------------
% PUTILS
%-----------------------------

cmd = sprintf('rm -f %s/putils.html',doc_dir);
system(cmd);

putils_fct = {'putils/mpistart.m'...
              'putils/mpirun.m'...
              'putils/mpiruns.m'};

m2html('mfiles',putils_fct,...
       'recursive','on',...
       'globalHypertextLinks','on',...
       'htmldir',doc_dir,...
       'source','on', ...
       'template',tmpl,...
       'indexFile','putils');

%--------------------------------

% Final clean-up in the stb_guide directory

cd(doc_dir);

% Fix links in top level files

old_str = '<a href="../sundialsTB.html">sundialsTB</a> &gt; <a href="../foo.html">foo</a>';
new_str = '<a href="sundialsTB.html">sundialsTB</a>';

cmd = sprintf('sed ''s$%s$%s$'' install_STB.html > tmp_file',old_str,new_str);
system(cmd);
system('mv tmp_file install_STB.html');

%cmd = sprintf('sed ''s$%s$%s$'' startup_STB.html > tmp_file',old_str,new_str);
%system(cmd);
%system('mv tmp_file startup_STB.html');


% Fix paths to images

system('sed ''s$../../../$../../$g'' install_STB.html > tmp_file');
system('mv tmp_file install_STB.html');

%system('sed ''s$../../../$../../$g'' startup_STB.html > tmp_file');
%system('mv tmp_file startup_STB.html');


% Remove generated files not needed

system('rm -f foo.html');
system('rm -f cvodes/cvodes.html');
system('rm -f cvodes/cvm/cvodes.html');
system('rm -f cvodes/examples_ser/cvodes.html');
system('rm -f cvodes/examples_par/cvodes.html');
system('rm -f idas/idas.html');
system('rm -f idas/cvm/idas.html');
system('rm -f idas/examples_ser/idas.html');
system('rm -f idas/examples_par/idas.html');
system('rm -f kinsol/kinsol.html');
system('rm -f kinsol/kim/kinsol.html');
system('rm -f kinsol/examples_ser/kinsol.html');
system('rm -f kinsol/examples_par/kinsol.html');
system('rm -f putils/putils.html');
system('rm -f nvector/nvector.html');

% Overwrite files for which we have nicer ones :-)

system('rm -f cvodes.html idas.html kinsol.html nvector.html putils.html');

system('cp ../html_files/cvodes_top.html cvodes.html');
system('cp ../html_files/idas_top.html idas.html');
system('cp ../html_files/kinsol_top.html kinsol.html');
system('cp ../html_files/nvector_top.html nvector.html');
system('cp ../html_files/putils_top.html putils.html');

% Add other needed files

system('cp ../html_files/sundialsTB.html .');
system('cp ../html_files/*.png .');

cd('..');
