
"""
   runtest - numerically tolerant test library.

   Author:
      Radovan Bast (lastname at kth.se).

   License:
      GNU Lesser General Public License.

   Documentation:
      http://runtest.readthedocs.org

   Source:
      https://github.com/rbast/runtest
"""

RUNTEST_VERSION = 'v0.1.3'

import re
import os
import sys
import subprocess
import shlex
import shutil
import string
from optparse import OptionParser


#------------------------------------------------------------------------------
class FilterKeywordError(Exception):
    pass


#------------------------------------------------------------------------------
class TestFailedError(Exception):
    pass


#------------------------------------------------------------------------------
class BadFilterError(Exception):
    pass


#------------------------------------------------------------------------------
class AcceptedError(Exception):
    pass


#------------------------------------------------------------------------------
class SubprocessError(Exception):
    pass


#------------------------------------------------------------------------------
class TestRun:

    #--------------------------------------------------------------------------
    def __init__(self, _file, argv):

        self.input_dir = input_dir = os.path.dirname(os.path.realpath(_file))

        options = self._parse_args(input_dir, argv)
        self.binary_dir = options.binary_dir
        self.work_dir   = options.work_dir
        self.verbose    = options.verbose
        self.skip_run   = options.skip_run

        if self.work_dir != self.input_dir:
            self._safe_copy(self.input_dir, self.work_dir)

        os.chdir(self.work_dir) # FIXME possibly problematic

    #--------------------------------------------------------------------------
    def execute(self,
                command,
                stdout_file_name='',
                accepted_errors=[]):
        """
        Runs the command.

        Raises:
            - AcceptedError
            - SubprocessError
        """
        if self.skip_run:
            return

        if sys.platform != "win32":
            command = shlex.split(command)

        process = subprocess.Popen(command,
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        for error in accepted_errors:
            if error in stderr:
                # we found an error that we expect/accept
                raise AcceptedError('found error which is expected/accepted: %s\n' % error)
        if process.returncode != 0:
            raise SubprocessError('ERROR: crash during %s\n%s' % (command, stderr))
        if stdout_file_name != '':
            f = open(stdout_file_name, 'w')
            f.write(stdout)
            f.close()

    #--------------------------------------------------------------------------
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

    #--------------------------------------------------------------------------
    def _parse_args(self, input_dir, argv):
        parser = OptionParser(description='runtest %s - Numerically tolerant test library.' % RUNTEST_VERSION)
        parser.add_option('--binary-dir',
                          '-b',
                          action='store',
                          default=input_dir,
                          help='directory containing the binary/launcher [default: %default]')
        parser.add_option('--work-dir',
                          '-w',
                          action='store',
                          default=input_dir,
                          help='working directory [default: %default]')
        parser.add_option('--verbose',
                          '-v',
                          action='store_true',
                          default=False,
                          help='give more verbose output upon test failure [default: %default]')
        parser.add_option('--skip-run',
                          '-s',
                          action='store_true',
                          default=False,
                          help='skip actual calculation(s) [default: %default]')
        (options, args) = parser.parse_args(args=argv[1:])

        if sys.platform == "win32":
            # on windows we flip possibly wrong slashes
            options.binary_dir = string.replace(options.binary_dir, '/', '\\')
            options.work_dir = string.replace(options.work_dir, '/', '\\')

        return options


#------------------------------------------------------------------------------
class _SingleFilter:

    def __init__(self, **kwargs):
        recognized_keywords = ['from_re',
                               'to_re',
                               're',
                               'from_string',
                               'to_string',
                               'string',
                               'ignore_below',
                               'ignore_sign',
                               'mask',
                               'num_lines',
                               'rel_tolerance',
                               'abs_tolerance']

        # check for unrecognized keywords
        for key in kwargs.keys():
            if key not in recognized_keywords:
                available_keywords = (', ').join(recognized_keywords)
                message = 'ERROR: keyword "%s" not recognized\n       ' % key
                message += 'available keywords: %s\n' % available_keywords
                raise FilterKeywordError(message)

        # check for incompatible keywords
        self._check_incompatible_keywords('from_re'      , 'from_string'  , kwargs)
        self._check_incompatible_keywords('to_re'        , 'to_string'    , kwargs)
        self._check_incompatible_keywords('to_string'    , 'num_lines'    , kwargs)
        self._check_incompatible_keywords('to_re'        , 'num_lines'    , kwargs)
        self._check_incompatible_keywords('string'       , 'from_string'  , kwargs)
        self._check_incompatible_keywords('string'       , 'to_string'    , kwargs)
        self._check_incompatible_keywords('string'       , 'from_re'      , kwargs)
        self._check_incompatible_keywords('string'       , 'to_re'        , kwargs)
        self._check_incompatible_keywords('string'       , 'num_lines'    , kwargs)
        self._check_incompatible_keywords('re'           , 'from_string'  , kwargs)
        self._check_incompatible_keywords('re'           , 'to_string'    , kwargs)
        self._check_incompatible_keywords('re'           , 'from_re'      , kwargs)
        self._check_incompatible_keywords('re'           , 'to_re'        , kwargs)
        self._check_incompatible_keywords('re'           , 'num_lines'    , kwargs)
        self._check_incompatible_keywords('rel_tolerance', 'abs_tolerance', kwargs)

        # now continue with keywords
        self.from_string = kwargs.get('from_string', '')
        self.to_string = kwargs.get('to_string', '')
        self.ignore_sign = kwargs.get('ignore_sign', False)
        self.ignore_below = kwargs.get('ignore_below', 1.0e-40)
        self.num_lines = kwargs.get('num_lines', 0)

        if 'rel_tolerance' in kwargs.keys():
            self.tolerance = kwargs.get('rel_tolerance')
            self.tolerance_is_relative = True
            self.tolerance_is_set = True
        elif 'abs_tolerance' in kwargs.keys():
            self.tolerance = kwargs.get('abs_tolerance')
            self.tolerance_is_relative = False
            self.tolerance_is_set = True
        else:
            self.tolerance_is_set = False

        self.mask = kwargs.get('mask', [])
        if self.mask == []:
            self.use_mask = False
        else:
            self.use_mask = True
            for i in self.mask:
                if i < 1:
                    raise FilterKeywordError('ERROR: mask starts counting from 1 (first word)\n')

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

        only_string = kwargs.get('string', '')
        if only_string != '':
            self.from_string = only_string
            self.num_lines = 1

        only_re = kwargs.get('re', '')
        if only_re != '':
            self.from_string = only_re
            self.num_lines = 1
            self.from_is_re = True

    def _check_incompatible_keywords(self, kw1, kw2, kwargs):
        if kw1 in kwargs.keys() and kw2 in kwargs.keys():
            raise FilterKeywordError('ERROR: incompatible keywords: "%s" and "%s"\n' % (kw1, kw2))


#------------------------------------------------------------------------------
class Filter:

    #--------------------------------------------------------------------------
    def __init__(self):
        self.filter_list = []

    #--------------------------------------------------------------------------
    def add(self, *args, **kwargs):
        """
        Adds filter task to list of filters.

        Raises:
            - FilterKeywordError
        """
        self.filter_list.append(_SingleFilter(*args, **kwargs))

    #--------------------------------------------------------------------------
    def check(self, work_dir, out_name, ref_name, verbose=False):
        """
        Compares output (work_dir/out_name) with reference (work_dir/ref_name)
        applying all filters tasks from the list of filters.

        Input:
            - work_dir -- working directory
            - out_name -- actual output file name
            - ref_name -- reference output file name
            - verbose  -- give verbose output upon failure

        Returns:
            - nothing

        Generates the following files in work_dir:

            - out_name.filtered  -- numbers extracted from output
            - out_name.reference -- numbers extracted from reference
            - out_name.diff      -- difference between the two above

        Raises:
            - TestFailedError
            - BadFilterError
        """

        out = open(out_name).readlines()
        log_out = open('%s.filtered' % out_name, 'w')

        ref = open(ref_name).readlines()
        log_ref = open('%s.reference' % out_name, 'w')

        log_diff = open('%s.diff' % out_name, 'w')

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

        pattern_int = re.compile('^-?[0-9]+$', re.VERBOSE)
        pattern_float = re.compile(numeric_const_pattern, re.VERBOSE)
        pattern_d = re.compile(r'[dD]')

        for f in self.filter_list:

            line_l_out = self._extract(f, out)
            line_l_ref = self._extract(f, ref)
            for (line_l, string_l, log, file_name) in \
                    [(line_l_out, out, log_out, out_name),
                     (line_l_ref, ref, log_ref, ref_name)]:
                if line_l == []:
                    if f.num_lines > 0:
                        r = '[%i lines from "%s"]' % (f.num_lines, f.from_string)
                    else:
                        r = '["%s" ... "%s"]' % (f.from_string, f.to_string)
                    message = 'ERROR: filter %s did not extract anything from file %s\n' % (r, file_name)
                    raise BadFilterError(message)
                for line in line_l:
                    log.write(string_l[line])

            f_l_out = []
            f_l_ref = []
            f_to_line_out = []
            f_to_line_ref = []
            for (line_l, string_l, f_l, f_to_line) in \
                    [(line_l_out, out, f_l_out, f_to_line_out),
                     (line_l_ref, ref, f_l_ref, f_to_line_ref)]:
                for line in line_l:
                    i = 0
                    for w in string_l[line].split():
                        i += 1
                        if (f.use_mask) and (i not in f.mask):
                            continue
                        # do not consider words like TzB1g
                        # otherwise we would extract 1 later
                        if re.match(r'^[0-9\.eEdD\+\-]*$', w):
                            is_integer = False
                            if len(pattern_float.findall(w)) > 0:
                                is_integer = (pattern_float.findall(w) == pattern_int.findall(w))
                            # apply floating point regex
                            for m in pattern_float.findall(w):
                                # substitute dD by e
                                m = pattern_d.sub('e', m)
                                if is_integer:
                                    f_l.append(int(m))
                                else:
                                    f_l.append(float(m))
                                f_to_line.append(line)

            if len(f_l_out) == len(f_l_ref):
                for i in range(len(f_l_out)):
                    r_out = f_l_out[i]
                    r_ref = f_l_ref[i]

                    if f.ignore_sign:
                        # if ignore sign take absolute values
                        r_out = abs(r_out)
                        r_ref = abs(r_ref)

                    is_integer_out = isinstance(r_out, int)
                    is_integer_ref = isinstance(r_ref, int)

                    if is_integer_out and is_integer_ref:
                        # we compare integers
                        if r_out != r_ref:
                            log_diff.write('line %i: %s' % (f_to_line_out[i] + 1, out[f_to_line_out[i]]))
                            log_diff.write('    found integer: %i, expected: %i\n\n' % (r_out, r_ref))
                    else:
                        # we compare floats
                        if not f.tolerance_is_set:
                            raise FilterKeywordError('ERROR: for floats you have to specify either rel_tolerance or abs_tolerance\n')
                        if abs(r_ref) > f.ignore_below:
                            # calculate relative error only for
                            # significant ('nonzero') numbers
                            error = r_out - r_ref
                            if f.tolerance_is_relative:
                                error /= r_ref
                            if abs(error) > f.tolerance:
                                log_diff.write('line %i: %s' % (f_to_line_out[i] + 1, out[f_to_line_out[i]]))
                                if f.tolerance_is_relative:
                                    log_diff.write('    rel error %7.4e > tolerance %7.4e\n\n' % (error, f.tolerance))
                                else:
                                    log_diff.write('    abs error %7.4e > tolerance %7.4e\n\n' % (error, f.tolerance))
            else:
                log_diff.write('extracted sizes do not match\n')

                log_diff.write('own gave %i numbers:\n' % len(f_l_out))
                last_line_printed = -1
                for i in range(len(f_l_out)):
                    if (f_to_line_out[i] != last_line_printed):
                        log_diff.write('    %s' % out[f_to_line_out[i]])
                        last_line_printed = f_to_line_out[i]

                log_diff.write('reference gave %i numbers:\n' % len(f_l_ref))
                last_line_printed = -1
                for i in range(len(f_l_ref)):
                    if (f_to_line_ref[i] != last_line_printed):
                        log_diff.write('    %s' % ref[f_to_line_ref[i]])
                        last_line_printed = f_to_line_ref[i]
        log_out.close()
        log_ref.close()
        log_diff.close()
        if os.path.getsize('%s.diff' % out_name) > 0:
            log_diff = open('%s.diff' % out_name, 'r')
            diff = ''
            for line in log_diff.readlines():
                diff += line
            log_diff.close()
            message = "ERROR: test %s failed\n" % out_name
            if verbose:
                message += diff
            raise TestFailedError(message)

    #--------------------------------------------------------------------------
    def _extract(self, f, out):
        """
        Input:
            - f -- filter task
            - out -- the output string to filter

        Returns:
            - line_l -- list of line numbers from which numbers will be extracted

        Raises:
            nothing
        """

        if f.from_is_re:
            ps = re.compile(r'.*%s' % f.from_string)
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
                if f.num_lines > 0:
                    for n in range(i, i + f.num_lines):
                        line_l.append(n)
                else:
                    for j in range(i, len(out)):
                        f.end_line_matches = False
                        if f.to_is_re:
                            f.end_line_matches = pe.match(out[j])
                        else:
                            f.end_line_matches = (f.to_string in out[j])
                        if f.end_line_matches:
                            for n in range(i, j + 1):
                                line_l.append(n)
                            break

        return line_l
