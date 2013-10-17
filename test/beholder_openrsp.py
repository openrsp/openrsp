
import os
import sys
from beholder import TestRun

class OpenRSPTestRun(TestRun):

    def __init__(self, _file, argv):
        TestRun.__init__(self, _file, argv)

    def execute(self, inp_files, mol_files, f=None, args=''):
        launcher = os.path.join(self.binary_dir, 'dalton')
        launcher += ' -noarch -nobackup'
        if args != '':
            launcher += ' %s' % args
        for inp in inp_files:
            inp_no_suffix = os.path.splitext(inp)[0]
            for mol in mol_files:
                mol_no_suffix = os.path.splitext(mol)[0]
                output_no_suffix = '%s_%s' % (inp_no_suffix, mol_no_suffix)
                sys.stdout.write('\nrunning test %s %s\n' % (inp_no_suffix, mol_no_suffix))
                TestRun.execute(self,
                                command = launcher + ' %s %s' % (inp, mol),
                                success_message = 'Total wall time used in DALTON',
                                success_file_name = '%s.out' % output_no_suffix,
                                stdout_file_name = '%s.stdout' % output_no_suffix,
                                stderr_means_crash = False)
                if f != None:
                    # f is a suffix-filter dictionary
                    for suffix in f:
                        out = '%s.%s' % (output_no_suffix, suffix)
                        f[suffix].check(self.work_dir, '%s' % out, 'result/%s' % out)
                sys.stdout.write('passed\n')
