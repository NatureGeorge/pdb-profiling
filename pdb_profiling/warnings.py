# @Created Date: 2020-10-26 01:00:07 pm
# @Filename: warnings.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-26 01:00:42 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import warnings

def custom_warn_format(msg, cat, filename, lineno, file=None, line=None):
    return f'{cat.__name__}: {msg} (from {filename}:{lineno})'


warnings.formatwarning = custom_warn_format

class MultiWrittenWarning(UserWarning):
    pass


class WithoutCifKeyWarning(UserWarning):
    pass

class PISAErrorWarning(UserWarning):
    pass


class MultipleConformersWarning(UserWarning):
    pass


class ConflictChainIDWarning(UserWarning):
    pass


class PossibleObsoletedUniProtWarning(UserWarning):
    pass


class PossibleObsoletedPDBEntryWarning(UserWarning):
    pass


class SkipAssemblyWarning(UserWarning):
    pass


class PeptideLinkingWarning(UserWarning):
    pass
