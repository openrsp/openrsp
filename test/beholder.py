
"""
    Beholder - numerically tolerant test library.

    Version X.Y

    Copyright (c) 2013, Radovan Bast (lastname at kth.se)
    Licensed under the GNU Lesser General Public License.

    For documentation see https://github.com/rbast/beholder
"""

import re
import os
import sys
import subprocess
import shlex
import shutil
from optparse import OptionParser

#-------------------------------------------------------------------------------
class TestRun():

    #---------------------------------------------------------------------------
    def __init__(self, _file, argv):
        self.input_dir = input_dir = os.path.dirname(os.path.realpath(_file))
        self.binary_dir, self.work_dir = self._parse_args(input_dir, argv)
        if self.work_dir != self.input_dir:
            self._safe_copy(self.input_dir, self.work_dir)
        os.chdir(self.work_dir)

    #---------------------------------------------------------------------------
    def execute(self,
                command,
                stdout_file_name='',
                accepted_errors=[]):
        try:
            process = subprocess.Popen(shlex.split(command),
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            for error in accepted_errors:
                if error in stderr:
                    # we found an error that we expect/accept
                    # in this case we do not mark test as red and exit with zero
                    sys.stderr.write('found error which is expected/accepted: %s\n' % error)
                    sys.exit(0)
            if process.returncode != 0:
                sys.stderr.write('ERROR: crash during %s\n' % command)
                sys.stderr.write(stderr)
                sys.exit(-1)
            if stdout_file_name != '':
                f = open(stdout_file_name, 'w')
                f.write(stdout)
                f.close()
        except OSError:
            sys.stderr.write('ERROR: could not execute command %s\n' % command)
            sys.stderr.write('       have you set the correct --binary-dir (or -b)?\n')
            sys.stderr.write('       try also --help\n')
            sys.exit(-1)

    #---------------------------------------------------------------------------
    def _safe_copy(self, root_src_dir, root_dst_dir, exclude_files=[]):
        for src_dir, dirs, files in os.walk(root_src_dir):
            dst_dir = src_dir.replace(root_src_dir, root_dst_dir)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            for f in files:
                if f not in exclude_files:
                    src_file = os.path.join(src_dir, f)
                    dst_file = os.path.join(dst_dir, f)
                    shutil.copy(src_file, dst_file)

    #---------------------------------------------------------------------------
    def _parse_args(self, input_dir, argv):
        parser = OptionParser(description='Beholder - numerically tolerant test library.')
        parser.add_option('--binary-dir',
                          '-b',
                          action='store',
                          default=input_dir,
                          help='directory containing the binary/launcher [default: %(default)s]')
        parser.add_option('--work-dir',
                          '-w',
                          action='store',
                          default=input_dir,
                          help='working directory [default: %(default)s]')
        (options, args) = parser.parse_args(args=argv[1:])
        return (options.binary_dir, options.work_dir)

#-------------------------------------------------------------------------------
class _SingleFilter():

    def __init__(self, **kwargs):

        self.from_string  = kwargs.get('from_string', '')
        self.to_string    = kwargs.get('to_string', '')
        self.tolerance    = kwargs.get('tolerance', 1.0e-5)
        self.ignore_sign  = kwargs.get('ignore_sign', False)
        self.ignore_below = kwargs.get('ignore_below', 1.0e-40)
        self.nr_lines     = kwargs.get('nr_lines', 0)
        self.mask         = kwargs.get('mask', [])
        self.use_mask     = False
        if self.mask != []:
            self.use_mask = True

        self.from_is_re = False
        from_re = kwargs.get('from_re', '')
        if from_re != '':
            self.from_string = from_re
            self.from_is_re = True

        self.to_is_re = False
        to_re = kwargs.get('to_re', '')
        if to_re != '':
            self.to_string = to_re
            self.to_is_re = True

#-------------------------------------------------------------------------------
class Filter():

    #---------------------------------------------------------------------------
    def __init__(self):
        self.filter_list = []

    #---------------------------------------------------------------------------
    def add(self, *args, **kwargs):
        self.filter_list.append(_SingleFilter(*args, **kwargs))

    #---------------------------------------------------------------------------
    def check(self, work_dir, out_name, ref_name):

        os.chdir(work_dir)

        log_diff = open('%s.diff'      % out_name, 'w')
        log_out  = open('%s.result'    % out_name, 'w')
        log_ref  = open('%s.reference' % out_name, 'w')

        out = open(out_name).readlines()

        if not os.path.exists(ref_name):
            sys.stderr.write('ERROR: reference output %s not found\n' % ref_name)
            sys.exit(-1)
        ref = open(ref_name).readlines()

        for f in self.filter_list:
            line_l     = self._extract(f, out, out_name, log_out)
            line_l_ref = self._extract(f, ref, ref_name, log_ref)

            numeric_const_pattern = r"""
            [-+]? # optional sign
            (?:
                (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                |
                (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
            )
            # followed by optional exponent part if desired
            (?: [EeDd] [+-]? \d+ ) ?
            """

            p  = re.compile(numeric_const_pattern, re.VERBOSE)
            pd = re.compile(r'[dD]')

            f_l           = []
            f_l_ref       = []
            f_to_line     = []
            f_to_line_ref = []
            for x in line_l:
                for line in x:
                    i = 0
                    for w in out[line].split():
                        i += 1
                        if (f.use_mask) and (i not in f.mask): break
                        # do not consider words like TzB1g
                        # otherwise we would extract 1 later
                        if re.match(r'^[0-9\.eEdD\+\-]*$', w):
                            # apply floating point regex
                            for m in p.findall(w):
                                 # substitute dD by e
                                 m = pd.sub('e', m)
                                 f_l.append(float(m))
                                 f_to_line.append(line)

            for x in line_l_ref:
                for line in x:
                    i = 0
                    for w in ref[line].split():
                        i += 1
                        if (f.use_mask) and (i not in f.mask): break
                        if re.match(r'^[0-9\.eEdD\+\-]*$', w):
                            for m in p.findall(w):
                                 m = pd.sub('e', m)
                                 f_l_ref.append(float(m))
                                 f_to_line_ref.append(line)

            if len(f_l) == len(f_l_ref):
                for i in range(len(f_l)):
                    r     = f_l[i]
                    r_ref = f_l_ref[i]

                    if f.ignore_sign:
                        # if ignore sign take absolute values
                        r     = abs(r)
                        r_ref = abs(r_ref)
                    error = r - r_ref

                    rel_error = error
                    if abs(r_ref) > f.ignore_below:
                        # calculate relative error only for significant ('nonzero') numbers
                        rel_error /= r_ref
                    if abs(rel_error) > f.tolerance:
                        log_diff.write('line %i: %s' % (f_to_line[i]+1, out[f_to_line[i]]))
                        log_diff.write('    rel error %7.4e > %7.4e\n\n' % (rel_error, f.tolerance))
            else:
                log_diff.write('compare floats sizes do not match\n')

                log_diff.write('own gave %i floats:\n' % len(f_l))
                last_line_printed = -1
                for i in range(len(f_l)):
                    if (f_to_line[i] != last_line_printed):
                       log_diff.write('    %s' % out[f_to_line[i]])
                       last_line_printed = f_to_line[i]

                log_diff.write('refence gave %i floats:\n' % len(f_l_ref))
                last_line_printed = -1
                for i in range(len(f_l_ref)):
                    if (f_to_line_ref[i] != last_line_printed):
                       log_diff.write('    %s' % ref[f_to_line_ref[i]])
                       last_line_printed = f_to_line_ref[i]
        log_diff.close()
        log_out.close()
        log_ref.close()
        if os.path.getsize('%s.diff' % out_name) != 0:
            sys.stderr.write('ERROR: test %s failed\n' % out_name)
            log_diff = open('%s.diff' % out_name, 'r')
            for line in log_diff.readlines():
                sys.stderr.write(line)
            log_diff.close()
            sys.exit(-1)

    #---------------------------------------------------------------------------
    def _extract(self, f, out, out_name, log):

        if f.from_is_re:
            ps = re.compile(r'.*%s' % f.from_string)
        if f.nr_lines == 0:
            if f.to_is_re:
                pe = re.compile(r'.*%s' % f.to_string)

        line_l = []
        for i in range(len(out)):
            start_line_matches = False
            if f.from_is_re:
                start_line_matches = ps.match(out[i])
            else:
                start_line_matches = (f.from_string in out[i])
            if start_line_matches:
                if f.nr_lines > 0:
                    line_l.append(range(i, i+f.nr_lines))
                else:
                    for j in range(i, len(out)):
                        f.end_line_matches = False
                        if f.to_is_re:
                            f.end_line_matches = pe.match(out[j])
                        else:
                            f.end_line_matches = (f.to_string in out[j])
                        if f.end_line_matches:
                            line_l.append(range(i, j+1))
                            break

        if line_l == []:
            if f.nr_lines > 0:
                sys.stderr.write('ERROR: filter [%i lines from "%s"] did not extract anything from file %s\n' % (f.nr_lines, f.from_string, out_name))
            else:
                sys.stderr.write('ERROR: filter ["%s" ... "%s"] did not extract anything from file %s\n' % (f.from_string, f.to_string, out_name))
            sys.exit(-1)

        for x in line_l:
            for line in x:
                log.write(out[line])

        return line_l
