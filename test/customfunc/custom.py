# @Created Date: 2020-06-01 07:04:11 pm
# @Filename: custom.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2020-06-01 07:05:14 pm
# @Copyright (c) 2020 MinghuiGroup, Soochow University
from re import compile as re_compile


pat_hgmd = re_compile(r"(.+\.[0-9]+):p\.([A-z]{1})([0-9]+)([A-z]{1})")
pat_cgi = re_compile(r"([A-z]{1})([0-9]+)([A-z]{1})")

def handle_hgmd_df(hgmd_df):
    '''
    usecols=['PROT']
    '''
    hgmd_df[['ftId', 'Ref', 'Pos', 'Alt']] = hgmd_df.apply(
        lambda x: pat_hgmd.search(x['PROT']).groups(), axis=1, result_type='expand')
    hgmd_df.drop(columns=['PROT'], inplace=True)
    hgmd_df.Pos = hgmd_df.Pos.astype(int)
    return hgmd_df.drop_duplicates().to_dict('records')


def handle_exac_df(exac_df):
    '''
    usecols=['ENST', 'aa_ref', 'aa_pos', 'aa_alt']
    '''
    columns = dict(
        zip(
            ('ENST', 'aa_ref', 'aa_pos', 'aa_alt'),
            ('ftId', 'Ref', 'Pos', 'Alt')
        )
    )
    exac_df.rename(columns=columns, inplace=True)
    exac_df.Pos = exac_df.Pos.astype(int)
    return exac_df.drop_duplicates().to_dict('records')


def handle_premptc(premptc_df):
    def split_df(df, colName, sep):
        return df.drop([colName], axis=1).join(df[colName].str.split(sep, expand=True).stack().reset_index(level=1, drop=True).rename(colName))
    '''
    usecols=['Transcript_ID', 'aaref', 'aapos', 'aaalt']
    '''
    columns = dict(
        zip(
            ('Transcript_ID', 'aaref', 'aapos', 'aaalt'),
            ('ftId', 'Ref', 'Pos', 'Alt')
        )
    )
    premptc_df.rename(columns=columns, inplace=True)
    premptc_df = split_df(premptc_df, 'ftId', ';').drop_duplicates()
    premptc_df.Pos = premptc_df.Pos.astype(int)
    return premptc_df.drop_duplicates().to_dict('records')


def handle_cgi_df(cgi_df):
    '''
    usecols=['gene', 'mutation_unp']
    '''
    cgi_df[['Ref', 'Pos', 'Alt']] = cgi_df.apply(
        lambda x: pat_cgi.search(x['mutation_unp']).groups(), axis=1, result_type='expand')
    cgi_df.rename(columns={'gene': 'ftId'}, inplace=True)
    cgi_df.Pos = cgi_df.Pos.astype(int)
    return cgi_df.drop_duplicates().to_dict('records')
