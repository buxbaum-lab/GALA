### NOTE This script includes functions adapted from the gnomAD QC repository (https://github.com/broadinstitute/gnomad_qc).
### Minor modifications have been made to suit the specific needs of this analysis.

import hail as hl

hl.init()

import pprint as pp

from typing import Any, Tuple, Union, Optional, List, Dict

#from .generic import *
import pickle
import copy
import uuid

def pLI(ht: hl.Table, obs: hl.expr.Int64Expression, exp: hl.expr.Float64Expression) -> hl.Table:
    """
    Compute the pLI score using the observed and expected variant counts.
    
    The output Table will include the following annotations:
        - pLI - Probability of loss-of-function intolerance; probability that transcript falls into 
            distribution of haploinsufficient genes
        - pNull - Probability that transcript falls into distribution of unconstrained genes
        - pRec - Probability that transcript falls into distribution of recessive genes

    :param ht: Input Table.
    :param obs: Expression for the number of observed variants on each gene or transcript in `ht`.
    :param exp: Expression for the number of expected variants on each gene or transcript in `ht`.
    :return: StructExpression for the pLI score.
    """
    last_pi = {'Null': 0, 'Rec': 0, 'LI': 0}
    pi = {'Null': 1 / 3, 'Rec': 1 / 3, 'LI': 1 / 3}
    expected_values = {'Null': 1, 'Rec': 0.463, 'LI': 0.089}
    ht = ht.annotate(_obs=obs, _exp=exp)

    while abs(pi['LI'] - last_pi['LI']) > 0.001:
        last_pi = copy.deepcopy(pi)
        ht = ht.annotate(
            **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
        ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
        ht = ht.annotate(**{k: ht[k] / ht.row_sum for k, v in pi.items()})
        pi = ht.aggregate({k: hl.agg.mean(ht[k]) for k in pi.keys()})

    ht = ht.annotate(
        **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
    ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
    return ht.select(**{f'p{k}': ht[k] / ht.row_sum for k, v in pi.items()})

def oe_confidence_interval(ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression,
                           prefix: str = 'oe', alpha: float = 0.05, select_only_ci_metrics: bool = True) -> hl.Table:
    """
    Determine the confidence interval around the observed:expected ratio.
    
    For a given pair of observed (`obs`) and expected (`exp`) values, the function computes the density of the Poisson distribution
    (performed using Hail's `dpois` module) with fixed k (`x` in `dpois` is set to the observed number of variants) over a range of
    lambda (`lamb` in `dpois`) values, which are given by the expected number of variants times a varying parameter ranging between
    0 and 2. The cumulative density function of the Poisson distribution density is computed and the value of the varying parameter
    is extracted at points corresponding to `alpha` (defaults to 5%) and 1-`alpha`(defaults to 95%) to indicate the lower and upper
    bounds of the confidence interval.
    
    Function will have following annotations in the output Table in addition to keys:
        - {prefix}_lower - the lower bound of confidence interval
        - {prefix}_upper - the upper bound of confidence interval

    :param ht: Input Table with the observed and expected variant counts for pLoF, missense, and synonymous variants.
    :param obs: Expression for the observed variant counts of pLoF, missense, or synonymous variants in `ht`.
    :param exp: Expression for the expected variant counts of pLoF, missense, or synonymous variants in `ht`.
    :param prefix: Prefix of upper and lower bounds, defaults to 'oe'.
    :param alpha: The significance level used to compute the confidence interval, defaults to 0.05.
    :param select_only_ci_metrics: Whether to return only upper and lower bounds instead of keeping all the annotations except `_exp`, defaults to True.
    :return: Table with the confidence interval lower and upper bounds.
    """
    ht = ht.annotate(_obs=obs, _exp=exp)
    oe_ht = ht.annotate(_range=hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000))
    oe_ht = oe_ht.annotate(_range_dpois=oe_ht._range.map(lambda x: hl.dpois(oe_ht._obs, oe_ht._exp * x)))

    oe_ht = oe_ht.transmute(_cumulative_dpois=hl.cumulative_sum(oe_ht._range_dpois))
    max_cumulative_dpois = oe_ht._cumulative_dpois[-1]
    oe_ht = oe_ht.transmute(_norm_dpois=oe_ht._cumulative_dpois.map(lambda x: x / max_cumulative_dpois))
    oe_ht = oe_ht.transmute(
        _lower_idx=hl.argmax(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x < alpha, x))),
        _upper_idx=hl.argmin(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x > 1 - alpha, x)))
    )
    oe_ht = oe_ht.transmute(**{
        f'{prefix}_lower': hl.if_else(oe_ht._obs > 0, oe_ht._range[oe_ht._lower_idx], 0),
        f'{prefix}_upper': oe_ht._range[oe_ht._upper_idx]
    })
    if select_only_ci_metrics:
        return oe_ht.select(f'{prefix}_lower', f'{prefix}_upper')
    else:
        return oe_ht.drop('_exp')
    

def add_rank(ht, field, ascending=True, total_genes=None, bins=10, defined_only=True):
    if total_genes is None:
        if defined_only:
            total_genes = ht.aggregate(hl.agg.count_where(hl.is_defined(ht[field])))
        else:
            total_genes = ht.count()
    rank_field = ht[field] if ascending else -ht[field]
    ht = ht.order_by(rank_field).add_index(f'{field}_rank')
    ht = ht.annotate(**{f'{field}_rank': hl.or_missing(
        hl.is_defined(ht[field]), ht[f'{field}_rank']
    )})
    return ht.annotate(**{
        f'{field}_rank': hl.int(ht[f'{field}_rank']),
        f'{field}_bin': hl.int(ht[f'{field}_rank'] * bins / total_genes),
        f'{field}_bin_6': hl.int(ht[f'{field}_rank'] * 6 / total_genes)
    })

ds = hl.import_table('gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz', force_bgz = True,
                    types={'exp_syn': hl.tfloat64, 'obs_syn': hl.tint64, 
                           'exp_lof': hl.tfloat64, 'obs_lof': hl.tint64, 
                           'exp_mis': hl.tfloat64, 'obs_mis': hl.tint64, 
                           'n_sites': hl.tint64, 'canonical': hl.tbool, 
                           'downsampling': hl.tint64, 'caf': hl.tfloat64})

ds.write('gs://google_bucket/gnomad.v2.1.1.lof_metrics.downsamplings.ht', overwrite = True)

ht = hl.read_table('gs://google_bucket/gnomad.v2.1.1.lof_metrics.downsamplings.ht')
ht.count() # 10199700

ht = ht.filter(ht.canonical, keep = True)
ht = ht.key_by('gene', 'transcript')
ht.count()

def calc_metrics(ht: hl.Table, obs: hl.expr.Int64Expression, exp: hl.expr.Float64Expression, pop: str = 'global',
                 sample_size: int = 125748, alpha: float = 0.05, select_only_ci_metrics: bool = True, 
                 ht_overwrite: bool = True, prefix: str = 'oe', 
                 path: str = 'gs://google_bucket/') -> hl.Table:
    
    """Calculates pLI and LOEUF score for given population at given sample size. Writes and exports .ht and .txt.bgz
    files with data."""
    
    tmp = ht.filter((ht.pop == pop) & (ht.downsampling == sample_size), keep = True)
    # calculate metrics
    pli_tmp = pLI(tmp, tmp.obs_lof2, tmp.exp_lof)
    loeuf_tmp = oe_confidence_interval(tmp, tmp.obs_lof2, tmp.exp_lof)
    # add binning
    pli_tmp = add_rank(pli_tmp, 'pLI')
    loeuf_tmp = add_rank(loeuf_tmp, 'oe_upper')
    # rename
    pli_tmp = pli_tmp.rename({'pLI':f'pLI_{pop}','pNull':f'pNULL_{pop}','pRec':f'pRec_{pop}',
                              'pLI_bin':f'pLI_bin_{pop}','pLI_bin_6':f'pLI_bin6_{pop}', 
                              'pLI_rank':f'pLI_rank_{pop}'})
    loeuf_tmp = loeuf_tmp.rename({'oe_lower':f'oe_lower_{pop}','oe_upper':f'LOEUF_{pop}',
                                  'oe_upper_bin':f'LOEUF_bin_{pop}','oe_upper_bin_6':f'LOEUF_bin6_{pop}',
                                  'oe_upper_rank':f'LOEUF_rank_{pop}'})   
    # add keys bc they got lost in the ranking
    pli_tmp = pli_tmp.key_by('gene', 'transcript')
    loeuf_tmp = loeuf_tmp.key_by('gene', 'transcript')
    # join!
    pli_tmp = pli_tmp.join(loeuf_tmp, how = 'outer')
    tmp = tmp.join(pli_tmp, how = 'outer')
    # write and export
    tmp.write(f'{path}constraint_metrics_{pop}_downsample_{sample_size}.ht', overwrite = ht_overwrite)
    tmp = hl.read_table(f'{path}constraint_metrics_{pop}_downsample_{sample_size}.ht').export(f'{path}constraint_metrics_{pop}_downsample_{sample_size}.txt.bgz')
    return tmp
    
ht = hl.read_table('gs://google_bucket/gnomad.v2.1.1.lof_metrics.downsamplings_canonical.ht')
ht.count()

# ALL of ALL size
glob = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'global', sample_size = 125748)
# global of NFE size
glob = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'global', sample_size = 56885)
# global of AMR size
glob = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'global', sample_size = 17296)
# global of AFR size
glob = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'global', sample_size = 8128)

# NFE of NFE size
nfe = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'nfe', sample_size = 56885)
# NFE of AMR size
nfe = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'nfe', sample_size = 17296)
# NFE of AFR size
nfe = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'nfe', sample_size = 8128)

# AMR of AMR size
amr = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'amr', sample_size = 17296)
# AMR of AFR size
amr = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'amr', sample_size = 8128)

# AFR of AFR size
afr = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'afr', sample_size = 8128)

# EAS of AFR size
eas = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'eas', sample_size = 8128)
# SAS of AFR size
sas = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'sas', sample_size = 8128)
# EAS of EAS size
eas = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'eas', sample_size = 9197)
# EAS of EAS size
sas = calc_metrics(ht, ht.obs_lof2, ht.exp_lof, pop = 'sas', sample_size = 15308)