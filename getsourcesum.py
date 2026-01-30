#!/usr/bin/env python3
# This py script extracts 'source' feature from annotation file.
import polars as pl
import sys
from dictdiffer import diff
from Bio import Entrez

Entrez.email = "w3const@g.nig.ac.jp"  # needs email address to use Entrez

def get_taxid_and_rank(scientific_name):
    handle = Entrez.esearch(db="taxonomy", term=scientific_name)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"] and len(record['IdList']) == 1:
        taxid = record["IdList"][0]
        fetch_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        fetch_record = Entrez.read(fetch_handle)
        fetch_handle.close()
        rank = fetch_record[0]["Rank"] if "Rank" in fetch_record[0] else ""
        return taxid, rank
    if record["IdList"] and len(record['IdList']) > 1:
        return 'more than one id found', 'more than one id found'
    if not record["IdList"]:
        return None, None

if __name__ == "__main__":
    args = sys.argv
    df = pl.read_csv(
        args[1],
        separator='\t',
        has_header=False,
        new_columns=['col1', 'col2', 'col3', 'col4', 'col5']
    ).with_row_index("row")

    rows_first_col = df.filter((pl.col('col1').is_not_null()) & (pl.col('col1') != '')).select('row').to_series().to_list()

    # find datatype
    dtypes = df.filter(pl.col('col2') == 'DATATYPE').select('col5').to_series().to_list()
    if len(dtypes) > 0:
        print(dtypes)
    else:
        print("No DATATYPE found")

    countsrc = 0
    for k, start in enumerate(rows_first_col):
        if k == len(rows_first_col) - 1:
            entry_df = df.slice(start, df.height - start)
        else:
            entry_df = df.slice(start, rows_first_col[k+1] - start)
        features_idx = entry_df.filter(pl.col('col2').is_not_null()).select('row').to_series().to_list()
        queryft_idx = entry_df.filter(pl.col('col2') == 'source').select('row').to_series().to_list()
        # print(features_idx, srces_idx)
        topology_idx = entry_df.filter(pl.col('col2') == 'TOPOLOGY').select('col4').to_series().to_list()
        if topology_idx and entry_df['col1'][0] == 'COMMON':
            print(entry_df['col1'][0], 'has TOPOLOGY', topology_idx[0])
        if not queryft_idx:
            print(entry_df['col1'][0], 'has No source feature')
            continue
        else:
            query_features_idx = [features_idx.index(i) for i in queryft_idx]        
        if not query_features_idx:
            print('No query-feature found')
        else:
            for v in query_features_idx:
                countsrc += 1
                if v + 1 == len(features_idx):
                    # (
                    # entry_df.filter((pl.col('row') >= features_idx[v]))
                    # .select(['col1', 'col2', 'col3', 'col4', 'col5'])
                    # .write_csv(sys.stdout, separator='\t', include_header=False)
                    # )
                    source_feat = (
                    entry_df.filter((pl.col('row') >= features_idx[v]))
                    .select(['col4', 'col5'])
                    )
                    dict_source_feat = dict(zip(source_feat['col4'].to_list(), source_feat['col5'].to_list()))
                else:
                    # (
                    # entry_df.filter((pl.col('row') >= features_idx[v]) & (pl.col('row') < features_idx[v+1]))
                    # .select(['col1', 'col2', 'col3', 'col4', 'col5'])
                    # .write_csv(sys.stdout, separator='\t', include_header=False)
                    # )
                    source_feat = (
                    entry_df.filter((pl.col('row') >= features_idx[v]) & (pl.col('row') < features_idx[v+1]))
                    .select(['col4', 'col5'])
                    )
                    dict_source_feat = dict(zip(source_feat['col4'].to_list(), source_feat['col5'].to_list()))
                    # TOPOLOGY if exists
                    if topology_idx:
                        dict_source_feat['TOPOLOGY'] = topology_idx[0]
                if countsrc == 1:
                    dict_firstsource_feat = dict_source_feat
                # print(countsrc)
                diffs = list(diff(dict_firstsource_feat, dict_source_feat))
                if diffs:
                    print('Differences with the 1st entry', entry_df['col1'][0], diffs)
    print('source of the 1st entry\n', dict_firstsource_feat)
    orgn = dict_firstsource_feat.get('organism')
    taxid, rank = get_taxid_and_rank(orgn)
    print("taxid:{} rank:{}".format(taxid, rank))
    print(countsrc, 'entries')