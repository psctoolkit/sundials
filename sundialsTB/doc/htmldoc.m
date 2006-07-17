% HTMLDOC - Creates html documentation for sundialsTB
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.6 $Date: 2006/03/26 21:34:51 $

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
              'cvodes/CVodeSetFSAOptions.m'...
              'cvodes/CVodeMalloc.m'...
              'cvodes/CVodeSensMalloc.m'...
              'cvodes/CVadjMalloc.m'...
              'cvodes/CVodeMallocB.m'...
              'cvodes/CVode.m'...
              'cvodes/CVodeB.m'...
              'cvodes/CVodeGetStats.m'...
              'cvodes/CVodeGetStatsB.m'...
              'cvodes/CVodeGet.m'...
              'cvodes/CVodeFree.m'...
              'cvodes/CVBandJacFn.m'...
              'cvodes/CVDenseJacFn.m'...
              'cvodes/CVGcommFn.m'...
              'cvodes/CVGlocalFn.m'...
              'cvodes/CVMonitorFn.m'...
              'cvodes/CVQuadRhsFn.m'...
              'cvodes/CVRhsFn.m'...
              'cvodes/CVRootFn.m'...
              'cvodes/CVSensRhsFn.m'...
              'cvodes/CVJacTimesVecFn.m'...
              'cvodes/CVPrecSetupFn.m'...
              'cvodes/CVPrecSolveFn.m'};

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
            'idas/IDASetFSAOptions.m'...
            'idas/IDAMalloc.m'...
            'idas/IDASensMalloc.m'...
            'idas/IDAadjMalloc.m'...
            'idas/IDAMallocB.m'...
            'idas/IDASolve.m'...
            'idas/IDASolveB.m'...
            'idas/IDAGetStats.m'...
            'idas/IDAGetStatsB.m'...
            'idas/IDAGet.m'...
            'idas/IDAFree.m'...
            'idas/IDABandJacFn.m'...
            'idas/IDADenseJacFn.m'...
            'idas/IDAGcommFn.m'...
            'idas/IDAGlocalFn.m'...
            'idas/IDAMonitorFn.m'...
            'idas/IDAQuadRhsFn.m'...
            'idas/IDAResFn.m'...
            'idas/IDARootFn.m'...
            'idas/IDASensResFn.m'...
            'idas/IDAJacTimesVecFn.m'...
            'idas/IDAPrecSetupFn.m'...
            'idas/IDAPrecSolveFn.m'};

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

kinsol_fct = {'kinsol/KINDenseJacFn.m'...
              'kinsol/KINBandJacFn.m'...
              'kinsol/KINFree.m'...
              'kinsol/KINGcommFn.m'...
              'kinsol/KINGetStats.m'...
              'kinsol/KINGlocalFn.m'...
              'kinsol/KINJacTimesVecFn.m'...
              'kinsol/KINMalloc.m'...
              'kinsol/KINPrecSetupFn.m'...
              'kinsol/KINPrecSolveFn.m'...
              'kinsol/KINSetOptions.m'...
              'kinsol/KINSol.m'...
              'kinsol/KINSysFn.m'};
 
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
