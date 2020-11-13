# @Created Date: 2020-11-13 02:35:34 pm
# @Filename: __init__.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-11-13 02:35:39 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
import orm
from unsync import unsync
from pdb_profiling.processors.database import SqliteDB


class UniProtDB(SqliteDB):

    def init_table_model(self):
        class VAR_SEQ(orm.Model):
            __tablename__ = 'VAR_SEQ'
            __metadata__ = self.metadata
            __database__ = self.database
            AltID = orm.String(max_length=100, primary_key=True, allow_blank=True)
            Entry = orm.String(max_length=50, primary_key=True)
            AltRange = orm.JSON(allow_null=True)
            AltInfo = orm.Text(allow_null=True, allow_blank=True)
            AltIso = orm.String(max_length=100, allow_blank=True, allow_null=True)
            AltLen = orm.JSON(allow_null=True)
            evidence = orm.Text(allow_null=True, allow_blank=True)

        self.VAR_SEQ = VAR_SEQ
