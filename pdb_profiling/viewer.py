# @Created Date: 2020-10-21 09:09:05 pm
# @Filename: viewer.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-10-21 09:09:10 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from pandas import isna, DataFrame
from pdb_profiling.utils import expand_interval
from pdb_profiling.processors.pdbe.record import PDB


class NGL(object):

    @staticmethod
    def get_related_auth_res(res_df, struct_asym_id, res_num_set):
        if len(res_num_set) > 0:
            subset = res_df[res_df.struct_asym_id.eq(struct_asym_id) & res_df.residue_number.isin(res_num_set)]
            return ('('+subset.author_residue_number.astype(str)+'^'+subset.author_insertion_code+subset.multiple_conformers.apply(lambda x: ' and %' if isna(x) else ' and %A')+')')
        else:
            subset = res_df[res_df.struct_asym_id.eq(struct_asym_id)]
            assert len(subset) == 1
            return ('('+subset.author_residue_number.astype(str)+'^'+subset.author_insertion_code+subset.multiple_conformers.apply(lambda x: ' and %' if isna(x) else ' and %A')+')').iloc[0]

    @staticmethod
    def get_type(x):
        if x in (
            'polypeptide(L)',
            'polypeptide(D)'):
            return 'protein'
        elif x in (
            'polydeoxyribonucleotide',
            'polyribonucleotide',
            'polydeoxyribonucleotide/polyribonucleotide hybrid'):
            return 'nucleic'
        else:
            return 'ligand'

    @classmethod
    def interface_sele_str_unit(cls, res_df, record, suffix):
        mol_type = cls.get_type(record[f'molecule_type{suffix}'])
        if mol_type != 'ligand':
            res = cls.get_related_auth_res(res_df, record[f'struct_asym_id{suffix}'], frozenset(expand_interval(record[f'interface_range{suffix}'])))
            sres = cls.get_related_auth_res(res_df, record[f'struct_asym_id{suffix}'], frozenset(expand_interval(record[f'surface_range{suffix}'])))
            res = ' or '.join(res)
            sres = ' or '.join(sres)
            return mol_type, f' and ({res})', f' and ({sres})'
        else:
            res = cls.get_related_auth_res(res_df, record[f'struct_asym_id{suffix}'], frozenset())
            sres = cls.get_related_auth_res(res_df, record[f'struct_asym_id{suffix}'], frozenset())
            return mol_type, f' and {res}', f' and {sres}'

    @classmethod
    def get_interface_sele_str(cls, record):
        res_df = PDB(record['pdb_id']).fetch_from_pdbe_api('api/pdb/entry/residue_listing/', PDB.to_dataframe).result()
        type_1, i_str_1, s_str_1 = cls.interface_sele_str_unit(res_df, record, '_1')
        type_2, i_str_2, s_str_2 = cls.interface_sele_str_unit(res_df, record, '_2')
        chain_id_1 = record['chain_id_1']
        chain_id_2 = record['chain_id_2']
        return ((f'{type_1} and :{chain_id_1}'+i_str_1, f'{type_2} and :{chain_id_2}'+i_str_2),
                (f'{type_1} and :{chain_id_1}'+s_str_1, f'{type_2} and :{chain_id_2}'+s_str_2))

    @classmethod
    def get_interface_view(cls, view, record, **kwargs):
        (i1, i2), (s1, s2) = cls.get_interface_sele_str(record)
        view.add_spacefill(selection=s1, opacity=kwargs.get('surface_opacity_1', 0.05), color=kwargs.get('surface_color_1', 'white'))
        view.add_spacefill(selection=s2, opacity=kwargs.get('surface_opacity_2', 0.05), color=kwargs.get('surface_color_2', 'white'))
        view.add_spacefill(selection=i1, opacity=kwargs.get('interface_opacity_1', 0.5), color=kwargs.get('interface_color_1', 'green'))
        view.add_spacefill(selection=i2, opacity=kwargs.get('interface_opacity_2', 0.5), color=kwargs.get('interface_color_2', 'red'))
        view.background = kwargs.get('background', '#F3F3F3')
        view._set_size(*kwargs.get('size', ('50%', '50%')))
        return view
