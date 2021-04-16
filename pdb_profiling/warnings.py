# @Created Date: 2020-10-26 01:00:07 pm
# @Filename: warnings.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-26 01:00:42 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import warnings
from rich.console import Console
from pathlib import Path

console = Console()
_warnings_showwarning = warnings.showwarning


def custom_warn_format(msg, category, filename, lineno, file=None, line=None):
    return f'[bold red]{category.__name__}[/bold red]: {msg} ({Path(filename)}:{lineno})'


warnings.formatwarning = custom_warn_format


def custom_showwarning(message, category, filename, lineno, file=None, line=None):
    if file is not None:
        _warnings_showwarning(message, category, filename, lineno, file, line)
    else:
        console.print(warnings.formatwarning(message, category, filename, lineno, line))


warnings.showwarning = custom_showwarning

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


class PossibleObsoletedUniProtIsoformWarning(UserWarning):
    pass


class ZeroSizeWarning(UserWarning):
    pass


class FileExistsWarning(UserWarning):
    pass


class InvalidFileContentWarning(UserWarning):
    pass


class WithoutRCSBClusterMembershipWarning(UserWarning):
    pass


class SequenceConflictWarning(UserWarning):
    pass


class PDBeKBResidueMappingErrorWarning(UserWarning):
    pass
