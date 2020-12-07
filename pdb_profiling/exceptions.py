# @Created Date: 2020-10-23 10:02:23 pm
# @Filename: exceptions.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-23 10:02:31 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University

class WithoutExpectedKeyError(Exception):
    pass


class PossibleConnectionError(Exception):
    pass


class PossibleObsoletedUniProtError(Exception):
    pass
