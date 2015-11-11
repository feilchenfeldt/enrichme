#!/usr/bin/env python
"""
TODO:
-- fix the gene summary mode so that is supports overlapping features
    (necessary to implement max_dist argument)
    DONE?
-- allow prior windowing of input data for better performance
-- ask magnus about applications


In top scores enrichment:
If data is on a grid: do not keep the values for each grid point,
but keep the top values and their distance.
Could this speed up things when shifting? But can be problematic
if real positions very clustered.

Speed up things:
Instead of searching features for each shift, make a big sparse matrix 
(or a dictionary) that says for each snp with which features it is associated.
Uemit suggest to use a bloomfilter! (fast dictionary)
 

-- Allow options to automatically run as many permutations as necessary 
to get something:
A) significant
B) significant above multiple testing


Make the output column names generic, i.e., independent
of the input (score, feature, category, n_features, ...)!?

"""
import os, sys
import pandas as pd
import numpy as np
import gc, logging
import multiprocessing as mp
from hs_vervet.tools import hs_pandas as hp
#import warnings
#warnings.simplefilter(action = "ignore", category = 'SettingWithCopyWarning')
#warnings.simplefilter("ignore")
eu = os.path.expanduser
jn = os.path.join
logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.DEBUG)


def get_sep(fn_fh):
    try:
        fn = fn_fh.name
    except AttributeError:
        fn = fn_fh
    ext =  os.path.splitext(fn)[-1]
    if ext == ".tsv":
        sep = "\t"
    elif ext == ".csv":
        sep = ","
    else:
        sep = None
        logging.warning('Automatically inferring file seperator of {}. '
                        'Consider renaming your file to *.tsv or *.csv for speed.'.format(fn))
    return sep

def init_rank_table(assoc):
    rt = pd.DataFrame({assoc.name:assoc.values,"rank":0,"out_of":0},index=assoc.index)
    rt.index.name = assoc.index.name
    return rt

def shift_rod(rod_df, rnd, mode = "grid"):
    """
    shift reference ordered data across the whole genome

    Input:
    rod_df ... pandas dataframe or series with mulitiindex (chrom, pos)

    modes ...
        'grid' ... just rotate the index of the rod data frame
                    this means that the positions stay the same only the
                    value for each position becomes different
                    Faster, but means that you only hit the same grid-point
                    this should make it conservative on large grids. Large
                    grids are problematic if the fraction of top windows 
                    considered becomes large 
        'continuous' ... add the random shift to each index value.
                         NOT IMPLEMENTED
    """
    if mode == "grid":
        new_start_i = int(len(rod_df)*rnd)
        rotate_data = np.concatenate((rod_df.iloc[new_start_i:].values,rod_df.iloc[:new_start_i].values))
        if  isinstance(rod_df,pd.core.series.Series):
            r = pd.Series(rotate_data,index=rod_df.index)
            return r
        elif isinstance(rod_df,pd.core.frame.DataFrame):
            r = pd.DataFrame(rotate_data,index=rod_df.index,columns=rod_df.columns)
            return r
    else:
        raise UserException("Only mode grid supported.")








#parallel support

def fun(f,q_in,q_out):
    while True:
        i,x = q_in.get()
        if i is None:
            break
        q_out.put((i,f(x)))

def parmap(f, X, nprocs):
    q_in   = mp.Queue(1)
    q_out  = mp.Queue()

    proc = [mp.Process(target=fun,args=(f,q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i,x in sorted(res)]


class CandidateEnrichment(object):
    """
    Test enrichment of a candidate feature list
    against a feature to category mapping.
    In a common example, candidate_features will be genes 
    and categories will gene ontology (GO) categories.
    
    Input:
    candidate_features ... list of feature names to be tested
    feature_to_category ... data frame with arbitrary index
    """

    def __init__(self, candidate_features, feature_to_category, value_name='candidate',
                 feature_name='feature', category_name='category', feature_df=None,  ncpus=1):

        self.candidate_features = np.unique(candidate_features)
        if len(self.candidate_features) < len(candidate_features):
            logging.warning("There are duplicates in candidate_features, "
                                                    "going to be removed...")

        self._bind_feature_to_category(feature_to_category, feature_name,  category_name)


        n_candidates = len(self.candidate_features)
        self.all_candidate_features = self.feature_to_category[self.feature_name].unique()
        self.candidate_features = np.intersect1d(self.candidate_features,
                                                 self.all_candidate_features, assume_unique=True)
        assert len(self.candidate_features) > 0, ("No candidate feature found in the "
                                   "first column of feature_to_category.")
        if len(self.candidate_features) < n_candidates:
            logging.warning("Only {} of the {} candidates are present "
                            "in feature_to_category mapping. "
                            "Only those are going to be used.".format(len(self.candidate_features), n_candidates))
        self.init_rank_table = self.initital_rank_table()
        self.rank_table = self.init_rank_table
        self.feature_df = feature_df
        self.ncpus = ncpus
        self.value_name = value_name


    def _bind_feature_to_category(self, feature_to_category, feature_name, category_name):
        assert feature_name in feature_to_category.columns
        assert category_name in feature_to_category.columns
        self.feature_to_category = feature_to_category
        self.feature_name = feature_name
        self.category_name = category_name



    def get_association(self, candidate_features):
        """
        Get series with number of candidate candidate_features associated with each
        category in feature_to_category.
        """
        try:
            assoc = self.feature_to_category.set_index(self.feature_name).ix[candidate_features].groupby(self.category_name).apply(len)
        except IndexError, e:
            raise e
        assoc.name = "n_" + self.feature_name
        assoc.index.name = self.category_name
        return assoc


    def initital_rank_table(self):
        logging.debug("get real assoc")
        real_assoc = self.get_association(self.candidate_features)
        logging.debug("init rank table")
        rt = init_rank_table(real_assoc)
        return rt

    def permuter(self, rank_table, n_permut):
        """
        Update the supplied rank table (rt) with the 
        """
        rt = rank_table.copy()
        for i in xrange(n_permut):
            permut_candidate_features = np.random.choice(self.all_candidate_features,
                                            size=len(self.candidate_features), replace=False)
            assoc = self.get_association(permut_candidate_features)
            rt["rank"] += (rt["n_candidate_features"] > assoc.reindex(rt.index).fillna(0))
            rt["out_of"] += 1
            gc.collect() #needed for low memory profile
        rt.sort('rank',ascending=False,inplace=True)
        return rt


    def get_permut_rank_table(self, n_permut):
        rti = self.init_rank_table
        if self.ncpus > 1:
            n_permut_proc = int(n_permut/self.ncpus)
            rts = parmap(lambda:self.permuter(rti,n_permut_proc), range(self.ncpus), self.ncpus)
            rt = reduce_mem(rts)
        else:
            rt = self.permuter(rti, n_permut)
        return rt

    def permute(self, n_permut):
        rt = self.get_permut_rank_table(n_permut)
        self.rank_table = reduce_mem([self.rank_table, rt])

    def get_pvals(self, pval_threshold=1, category_to_description=None):
        return rank_to_pval(self.rank_table, pval_threshold=1, category_to_description=None)

    def create_info(self):
        """
        Creates an attribute self.summary_per_feature.

        For CandidateEnrichment or TopScoresEnrichment 
        this data frame contains a boolean variable for
        each gene that is True if the gene is in the caniditate set.
        For SummaryEnrichment this attribute gives the summary 
        for each gene.
        """
        assert self.feature_df is not None, "Cannot create info without feature_df"
        summary_per_feature = pd.DataFrame({self.value_name:self.feature_df[self.feature_name].\
                                                                apply(lambda x: x in self.candidate_features).values},
                                                                        index=self.feature_df[self.feature_name].values)
        summary_per_feature[self.value_name] = summary_per_feature[self.value_name].astype(int)
        summary_per_feature.index.name = self.feature_name
        summary_per_feature['zscore'] = (summary_per_feature[self.value_name]-summary_per_feature[self.value_name].mean())\
                                        /summary_per_feature[self.value_name].std(ddof=0)
        summary_per_feature.sort('zscore', ascending=False, inplace=True)
        self.summary_per_feature = summary_per_feature



#    def get_candidate_location(self):
#        """
#        for a list of gene ids,
#        get a data frame with their.
#        position
#        """
#        try:
#            gi = self.feature_df[self.feature_df[self.feature_name].apply(lambda x: x in self.features)]
#        except Exception, e:
#            raise e
#        return gi



class SummaryEnrichment(CandidateEnrichment):
    """
    summary ... name of the function of the groupby object
                         to apply to the data (e.g. 'mean', 'max',...)
    max_dist ... is not implemented yet!!!


    two scenarios: 
    -- summary across features + summary across categories
    -- summary across features -> take top features 
        (can be modelled by more complicated summary function, eg, lambda x: mean(x)>thresh)? 
         -> count features per cat (summary with sum)
    """
    def __init__(self, value_s, feature_df,  feature_to_category, feature_name='feature', category_name='category',
                    feature_summary=None, feature_summary_fun=None,
                 category_summary=None, category_summary_fun=None, min_features_per_cat=2,
                                                        max_dist=0, chrom_len=None, ncpus=1):
        #todo: logging for AssertionErrors, explain what the problem is
        hp.check_rod(value_s) 
        hp.check_feature_df(feature_df)
        assert feature_name in feature_df.columns
        assert (feature_summary is not None) != (feature_summary_fun is not None), \
                                "Specify either feature_summary OR feature_summary_fun."
        assert (category_summary is not None) != (category_summary_fun is not None), \
                                "Specify either category_summary OR category_summary_fun."
        if chrom_len is not None:
            self.flat_index = True
        else:
            self.flat_index = False

        if self.flat_index:
            self.value_s = hp.rod_to_1d(value_s, chrom_len)
            self.feature_df = hp.feature_df_to_1d(feature_df, chrom_len)
        else:
            self.value_s = value_s.copy()
            self.feature_df = feature_df.copy()

        self.value_name = self.value_s.name

        self._bind_feature_to_category(feature_to_category, feature_name,  category_name)

        self.prune_feature_to_category(min_features_per_cat)

        self.feature_summary = feature_summary
        self.category_summary = category_summary
        self.feature_summary_fun = feature_summary_fun
        self.category_summary_fun = category_summary_fun
        self.max_dist = max_dist

        #logging.debug("I am before.")
        #import pdb
        #logging.debug("I am here.")
        #pdb.set_trace()
        #logging.debug("I am after.")

        self.init_rank_table = self.initital_rank_table()
        self.rank_table = self.init_rank_table
        self.ncpus = ncpus


    def prune_feature_to_category(self, min_features_per_cat):
        self.feature_to_category = self.feature_to_category.drop_duplicates()
        ftc_features = self.feature_to_category[self.feature_name].unique()
        features = self.feature_df[self.feature_name].unique()
        not_in_ftc =  np.setdiff1d(features, ftc_features, assume_unique=True)
        not_in_feature_df = np.setdiff1d(ftc_features,features, assume_unique=True)
        if len(not_in_ftc)>0:
            logging.warning("{} features from the features file are not "
                            "in the feature_to_category mapping.".format(len(not_in_ftc)))
        if len(not_in_feature_df)>0:
            logging.warning("{} features from the feature_to_category mapping are not "
                            "in the features file. Removing them.".format(len(not_in_feature_df)))
            self.feature_to_category = self.feature_to_category.set_index(self.feature_name).\
                                                        drop(pd.Index(not_in_feature_df)).reset_index()
        if min_features_per_cat > 0:
            n_cats = len(self.feature_to_category[self.category_name].unique())
            logging.info("Removing categories for which there are less than {} "
                            "features in the features file.".format(min_features_per_cat))
            self.feature_to_category = self.feature_to_category.groupby(self.category_name).\
                                                           filter(lambda x: len(x) >= min_features_per_cat)
            n_removed = n_cats - len(self.feature_to_category[self.category_name].unique())
            logging.info("{} categories removed.".format(n_removed))
        self.feature_to_category.reset_index(inplace=True)


    def initital_rank_table(self):
        real_assoc = self.get_association(self.value_s)
        rt = init_rank_table(real_assoc)
        return rt


    def get_association(self, value_s):
        summary_per_feature = self.get_summary_per_feature(value_s)
        assoc = self.get_summary_per_category(summary_per_feature)
        assoc.name = self.value_name + '_' + 'summary'# self.feature_name + self.feature_summary + '_' 
        assoc.index.name = self.category_name
        return assoc

    def get_summary_per_feature(self,value_s):
        if self.flat_index:
            values_per_feature = hp.data_per_feature_FI(value_s, self.feature_df,
                                                    feature_name=self.feature_name)
        else:
            values_per_feature = hp.data_per_feature(value_s, self.feature_df,
                                                    feature_name=self.feature_name, max_dist=self.max_dist)
        if self.feature_summary is not None:
            groups = values_per_feature.groupby(self.feature_name)
            return getattr(groups, self.feature_summary)()
        elif self.feature_summary_fun is not None:
            return self.feature_summary_fun(values_per_feature)
        #summary_per_feature = hp.apply_to_feature(values_per_feature,
        #                                          groupby_func_name=self.feature_summary,
        #                                         function=self.feature_summary_fun)
        #return summary_per_feature


    def get_summary_per_category(self,value_per_feature):
        """
        Calculates summary (e.g. mean) of values for the features
        in each of the given categories.

        Returns:
        series
        """
        value_to_category = self.feature_to_category.copy()
        values_per_feature_to_category = value_per_feature.ix[value_to_category[self.feature_name].values].values
        del value_to_category[self.feature_name]
        value_to_category[self.value_name] =  values_per_feature_to_category
        if self.category_summary is not None:
            groups = value_to_category.groupby(self.category_name)
            return getattr(groups, self.category_summary)()[self.value_name]
        elif self.category_summary_fun is not None:
            return self.category_summary_fun(value_to_category)[self.value_name]
        #summary_per_category = getattr(value_to_category.groupby('category'), self.category_summary)()
        #return summary_per_category['value']



    def permuter(self, rank_table, n_permut):
        rank_table = rank_table.copy()
        for rnd in np.random.rand(n_permut):
            s = shift_rod(self.value_s, rnd)
            assoc = self.get_association(s)
            rank_table["rank"] += (rank_table[assoc.name] > \
                                            assoc.reindex(rank_table.index).fillna(0))
            rank_table["out_of"] += 1
            gc.collect() #needed for low memory profile
        rank_table.sort('rank',ascending=False,inplace=True)
        return rank_table

    def create_info(self):
        summary_per_feature = self.get_summary_per_feature(self.value_s)
        #this is contained in the pvalue df
        #summary_per_category = self.get_summary_per_category(summary_per_feature)
        summary_per_feature['zscore'] = (summary_per_feature[self.value_name]-summary_per_feature[self.value_name].mean())\
                                            /summary_per_feature[self.value_name].std(ddof=0)
        summary_per_feature.sort('zscore', ascending=False, inplace=True)
        #self.summary_per_category = summary_per_category
        self.summary_per_feature = summary_per_feature



class TopScoresEnrichment(SummaryEnrichment):
    """
    Test for enrichment by shifting the
    input reference ordered scores.

    First, we search for features
    close to the top values of the reference 
    ordered input scores, and get associations
    between features and categories.

    Then, we repeat the following step many times:
    Randomly shift the input scores,
    against the (chrom, pos) index and
    repeat the first step for the shifted data.

    Repeating the second step many times,
    gives a null distribution of
    associations between features and categories,
    against which we compare the real assocations
    from step one.
    """
    def __init__(self,value_s, feature_df,  feature_to_category, top_type, top,
                                 feature_name='feature', category_name='category',
                                 min_features_per_cat=2, ascending=False, max_dist=0, ncpus=1):

        top_types = ['count','threshold','quantile']
        assert feature_name in feature_df.columns
        assert top_type in top_types, "top_type must be one of {}".format(top_types)
        self.value_s = value_s.copy()
        #the following prevents get_peaks from throwing unspecific error
        self.feature_df = feature_df.drop_duplicates(subset=feature_name)
        self._bind_feature_to_category(feature_to_category, feature_name,  category_name)
        self.value_name = self.value_s.name

        self.prune_feature_to_category(min_features_per_cat)

        if top_type == 'count':
            top_n = top
        elif top_type == 'quantile':
            top_n = int(len(self.value_s)*top)
        elif top_type == 'threshold':
            value_s_s = value_s.sort(ascending=ascending, inplace=False)
            if ascending:
                top_n = np.argmax(value_s_s.values>top)
            else:
                top_n = np.argmax(value_s_s.values<top)
            del value_s_s
        self.top_n = top_n
        self.ascending = ascending
        self.max_dist = max_dist
        self.init_rank_table = self.initital_rank_table()
        self.rank_table = self.init_rank_table
        self.ncpus = ncpus


    def get_association(self, value_s):
        value_s = value_s.sort(ascending=self.ascending, inplace=False)
        top_s = value_s.iloc[:self.top_n]
        del value_s
        candidate_features = hp.get_features(top_s, self.feature_df, feature_name=self.feature_name,
                                                                             max_dist=self.max_dist)
        assoc = CandidateEnrichment.get_association(self, candidate_features)
        #assoc = self.feature_to_category.set_index(self.feature_name).ix[cand_genes].groupby(self.category_name).apply(len)
        #assoc.name = "n_" + self.feature_name
        #assoc.index.name = self.category_name
        return assoc

    
    def get_peak_info(self,top_s,peaks_per_gene):
        peak_height_name = top_s.name
        gene_list_peak_pos = peaks_per_gene.reset_index([0])[self.feature_name].groupby(lambda x: x).apply(list)
        gene_list_peak_pos.name = "genes"
        gene_list_peak_pos.index = pd.MultiIndex.from_tuples(gene_list_peak_pos.index)
        peak_info = pd.concat([top_s,gene_list_peak_pos],axis=1)
        peak_info.sort(peak_height_name,ascending=False,inplace=True)
        peak_info.index.names = ["chrom","pos"]
        return peak_info

    def create_info(self):
        """
        Creates several info dataframe for the input data.
        self.peaks_per_feature .... top peaks per feature for features that
                                are hit by top peaks
                            
        self.top_peaks ... top peaks with the features that are close
        """
        value_s = self.value_s.sort(ascending=self.ascending, inplace=False)
        top_s = value_s.iloc[:self.top_n]
        self.candidate_features = hp.get_features(top_s, self.feature_df, feature_name=self.feature_name,
                                                                             max_dist=self.max_dist)
        sub_feature_df = self.feature_df.reset_index().set_index(self.feature_name).ix[self.candidate_features]\
                                                        .reset_index().set_index(['chrom','start'])
        self.peaks_per_feature = get_peaks(sub_feature_df, top_s, self.max_dist, feature_name=self.feature_name)
        features_sort_by_max = self.peaks_per_feature['peak_height'].groupby(lambda i:i[0]).max()
        features_sort_by_max.sort(ascending=False,inplace=True)
        self.peaks_per_feature = self.peaks_per_feature.ix[features_sort_by_max.index]
        self.top_peaks = self.get_peak_info(top_s, self.peaks_per_feature)
        #super(SummaryEnrichment, self).create_info()
        CandidateEnrichment.create_info(self)


def get_peaks(sub_gene_df,top_s,max_dist,feature_name):
    """
    For each gene in gene_info get the
    peaks within max_dist in top_s. This 
    is basically reverse engineering to get
    the peak info for each gene that was found 
    to be associated with a peak. 
    The reason for reverse engeneering rather than 
    storing this information when searching for the genes
    for each peak is that we want to use precisely the same
    function to search the genes for the real data and for the 
    permutations.


    Input:
    gene_info ... data frame with index ('chrom','start')
                and columns 'gene_id' and 'end'
    top_s ... series of peak positions with index (chrom, pos)
                and values peak height
    max_dist ... maximum distance between gene and peak
    """
    gene_info = sub_gene_df
    def get_dist(df,gene_pos):
        """
        calculate distance
        """
        s = pd.Series(df.index.droplevel(0).values - gene_pos.ix[df.index[0][0]],
                                                  index=df.index.droplevel(0).values)
        return s
    tot_gene_peaks_df = pd.DataFrame()
    if not top_s.index.is_monotonic:
        top_s = top_s.sortlevel([0,1])
    if not gene_info.index.is_monotonic:
        gene_info = gene_info.sort_index()
    for chrom in gene_info.index.droplevel(1).unique():
        loc_top_s = top_s.ix[chrom]
        start = np.searchsorted(loc_top_s.index.values+max_dist,gene_info.ix[chrom].index.values)
        end = np.searchsorted(loc_top_s.index.values-max_dist,gene_info.ix[chrom]["end"].values)
        x = pd.concat([loc_top_s.iloc[st:ed] for st,ed in zip(start,end)],
                          keys=gene_info.ix[chrom][feature_name].values)
        x.name = "peak_height"


        dist_start = x.groupby(lambda i: i[0]).\
                    apply(lambda df: get_dist(df,
                                              gene_info.ix[chrom].reset_index().set_index(feature_name)["start"]))
        dist_start.name = "dist_start"
        dist_end = x.groupby(lambda i: i[0]).\
                    apply(lambda df: get_dist(df,
                                              gene_info.ix[chrom].set_index(feature_name)["end"]))
        dist_end.name = "dist_end"
        gene_peaks_df = pd.concat([x,dist_start,dist_end],axis=1)
        gene_peaks_df.index = pd.MultiIndex.from_arrays([gene_peaks_df.index.droplevel(1),
                                         [chrom]*len(x),
                                         gene_peaks_df.index.droplevel(0)])
        tot_gene_peaks_df = pd.concat([tot_gene_peaks_df, gene_peaks_df],axis=0)
            

    tot_gene_peaks_df.index.names = [feature_name,"chrom","peak_pos"]
    return tot_gene_peaks_df


def get_peak_info(top_s,peaks_per_gene):
    peak_height_name = top_s.name
    gene_list_peak_pos = peaks_per_gene.reset_index([0])["gene_id"].groupby(lambda x: x).apply(list)
    gene_list_peak_pos.name = "genes"
    gene_list_peak_pos.index = pd.MultiIndex.from_tuples(gene_list_peak_pos.index)
    peak_info = pd.concat([top_s,gene_list_peak_pos],axis=1)
    peak_info.sort(peak_height_name,ascending=False,inplace=True)
    peak_info.index.names = ["chrom","pos"]
    return peak_info


def get_p_val(rank_table):
    """
    Input:
    
    """
    r =  1-rank_table["rank"]*1./(rank_table["out_of"]+1)
    r.sort()
    r.name = "p_value"
    return r




def open_reduce_fns(fns):
    permut_fhs = []
    for fn in fns:
        try:
            if os.stat(fn).st_size>0:
                permut_fhs.append(open(fn))
            else:
                logging.warning( "File seems to be empty. Skipping it: {}.".format(fn))
        except Exception, e:
            logging.warning("Can't open file. Skipping it: {}.".format(fn))
            logging.warning(str(e))
    return permut_fhs

def reduce_fhs(permut_fhs):
    tot_rank = pd.read_csv(permut_fhs[0], index_col=0, sep=get_sep(permut_fhs[0])).dropna()
    tot_rank["index"] = tot_rank.index
    tot_rank.drop_duplicates(subset="index",inplace=True)
    del tot_rank["index"]
    for fh in permut_fhs[1:]:
        rank_table = pd.read_csv(fh,index_col=0, sep=get_sep(fh)).dropna()
        rank_table["index"] = rank_table.index
        rank_table.drop_duplicates(subset="index",inplace=True)
        try:
            tot_rank["rank"] = tot_rank["rank"].add(rank_table["rank"],fill_value=0)
            tot_rank["out_of"] = tot_rank["out_of"].add(rank_table["out_of"],fill_value=0)
        except Exception, e:
            raise e
    return tot_rank

def reduce_mem(rank_tables):
    tot_rank = rank_tables[0].dropna()
    tot_rank["index"] = tot_rank.index
    tot_rank.drop_duplicates(subset="index",inplace=True)
    del tot_rank["index"]
    for rank_table in rank_tables[1:]:
        rank_table.dropna(inplace=True)
        rank_table["index"] = rank_table.index
        rank_table.drop_duplicates(subset="index",inplace=True)
        try:
            tot_rank["rank"] = tot_rank["rank"].add(rank_table["rank"],fill_value=0)
            tot_rank["out_of"] = tot_rank["out_of"].add(rank_table["out_of"],fill_value=0)
        except Exception, e:
            raise e
    return tot_rank

def rank_to_pval(rank_table, pval_threshold=1, category_to_description=None):
    try:
        rank_table.drop('p_value', axis=1, inplace=True)
    except ValueError:
        pass
    p_vals = get_p_val(rank_table)
    p_val_df = rank_table.join(p_vals)
    p_val_df.sort('p_value', inplace=True)
    p_val_df["benjamini_hochberg"] = p_val_df["p_value"] * \
                                    len(feature_to_category.iloc[:,1].unique())*1. /\
                                    np.arange(1,len(p_val_df)+1)

    p_val_df =  p_val_df[p_val_df["p_value"]<=pval_threshold]
    if category_to_description is not None:
        #ctd_cols = parse_cols(reduce_args.category_to_description_cols)
        #category_to_description = pd.read_csv(reduce_args.category_to_description,
        #                                 usecols=ctd_cols, sep='\t')
        category_to_description = category_to_description.\
                                    set_index(category_to_description.columns[0],drop=True)
        try:
            p_val_df.drop(category_to_description.columns[0], axis=1, inplace=True)
        except ValueError:
            pass
        p_val_df = p_val_df.join(category_to_description)
    #p_val_df.sort('p_value', inplace=True)
    return p_val_df

def save_info(enrich, name):
    enrich.create_info()
    enrich.summary_per_feature.to_csv(name + ".summary_per_feature.tsv", sep='\t')
    try:
        enrich.summary_per_category.to_csv(name + ".summary_per_category.tsv", sep='\t')
    except AttributeError:
        pass
    try:
        enrich.peaks_per_feature.to_csv(name + ".peaks_per_feature.tsv", sep='\t')
        enrich.top_peaks.to_csv(name + ".top_peaks.tsv", sep='\t')
    except AttributeError:
        pass

if __name__ == "__main__":
    import argparse, time
    #import pdb



    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description=
                                      "Test enrichment of features (e.g. genes) in certain categories.",
                                                                             add_help=False)
    parser.add_argument("-R",'--run_type', choices=['Permute','Reduce'])
    parser.add_argument('-N','--name', help='Base name for all the output files. ')
    parser.add_argument("--feature_to_category", type=argparse.FileType('r'),
                                            help="Filename for tsv that links features to"
                                                    " categories. E.g. go associations. ")
    parser.add_argument('--feature_to_category_cols',nargs=2, default = ['feature', 'category'],
                             help="Column labels or positions (0-indexed) of 'feature' and 'category'"
                                                        "In the feature_to_category tsv."
                                                        "Expects 2 integers or strings.")

    parser.add_argument('--logging_level','-l',
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                                default='INFO', help='Minimun level of logging.')
    parser.add_argument("--help",'-h', action='store_true',
                                                     help="Print help message and exit.")


    #Runtype parsers
    void_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    runtypeparsers = void_parser.add_subparsers(dest='run_type')

    permuteparser = runtypeparsers.add_parser('Permute', description='Make random permutations or shifts of input '
                                                              ' and save the ranks of the real assocations compared to permutations.')
    permuteparser.add_argument('-M','--mode', choices=['Candidate', 'Summary', 'TopScores'], help='Enrichment method.')
    permuteparser.add_argument('--n_permut', type=int, default=0, help='Number of permutations or shifts.')
                                                                    #'Use this option for multiple runs '
                                                                     #                'that are reduced in mode --M Reduce '
                                                                      #               'afterwards.')
    permuteparser.add_argument('--noinfo',action='store_true', help='Do not output additional info files such as '
                                                                     'summary values per gene and do not output pvalues. '
                                                                     'Use this option if you are running multiple permute '
                                                                     'job that are reduced later.')
    permuteparser.add_argument("--ncpus", '-nct',
                            type=int, default=1, help='Number of cpus to use in parallel.')

    #permuteparser.add_argument('--reduce', action='store_true', help='Reduce on the fly. Use this option if you '
    #                                                                 'just want to do a single job.')

    reduceparser = runtypeparsers.add_parser('Reduce', description='Reduces the results from multiple Permute-runs '
                                                              ' to a single table of enrichment pvalues.')
    reduceparser.add_argument('--permuts', nargs='*',
                                           help='Filenames of permutation results of individual runs. '
                                                ' (specified with --name in the permutation runs).')
    
    reduceparser.add_argument('--category_to_description',
                                 type=argparse.FileType('r'), help='Tsv with category to category description mapping.')
    reduceparser.add_argument('--category_to_description_cols', nargs=2, default = ['category', 'description'],
                                help="Column labels or positions (0-indexed) of 'category' and 'description' "
                                                        "In the category_to_description tsv. "
                                                        "Expects 2 integers or strings.")
    reduceparser.add_argument('--remove_input', action='store_true',
                                help="Delete the input files given by --permuts .")
    #reduceparser.add_argument('--pval_threshold', type=float, default=1,
    #                             help='Do not report categories with p_value>pval_threshold. '
    #                                  '(Of course, they are still used to calculate the Benjamini-Hochberg '
    #                                  ' multiple testing correction.)')

    #Mode parsers
    void_parser1 = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    modeparsers = void_parser1.add_subparsers(dest='mode')

    candidateparser = modeparsers.add_parser('Candidate', description='Test for enrichment of candidate features in feature categories.')
    candidateparser.add_argument('--candidates', type=argparse.FileType('r'),
                                                         help='Sequence of candidate genes.')

    
    summaryparser = modeparsers.add_parser('Summary', description='Test for enrichment by summarising values across features and '
                                                           'features across genes')
    topscoresparser = modeparsers.add_parser('TopScores', description='Test for enrichment of features that are close to top scores '
                                                               'in feature categories.')


    for p in [summaryparser, topscoresparser]:
        p.add_argument("--rod",type=argparse.FileType('r'),
                               help="Input reference ordered data tsv containing the values per position.")
        p.add_argument("--rod_cols", type=str, nargs=3, default = ['chrom','pos','value'],
                                            help="Column labels or positions (0-indexed) "
                                                        "of [chrom pos value]. "
                                                        "Expects 3 integers or strings.")

        p.add_argument("--features", type=argparse.FileType('r'),
                                        help="Filename of the feature info tsv "
                                             " with columns 'chrom' 'start' 'end' 'feature'.")
        p.add_argument('--feature_cols',nargs=4, default = ['chrom', 'start', 'end', 'feature'],
                             help="Column labels or positions (0-indexed) of ['chrom' 'start' 'end' 'feature']. "
                                                        "Expects 4 integers or  strings. "
                                                        "E.g., if the file is a bed with the "
                                                        "feature name in the 4th column, use --feature_cols 0 1 2 3")
        p.add_argument('--min_features_per_cat', type=int, default=2, help="Minimum number of features found "
                                                                           "in the features file in order not to exclude "
                                                                           "a category from the testing. If categories are "
                                                                           "excluded beforehand you can speed up calculations "
                                                                           "by setting this to 0.")
        #p.add_argument("--peaks_per_gene", type=argparse.FileType('w'),
        #                           help="File to write peak info for each gene.")
        #p.add_argument("--top_peaks", type=argparse.FileType('w'),
        #                                          help="File to write top peak info.")


    summaryparser.add_argument('--feature_summary', default='max', choices = ['mean', 'median', 'min', 'max', 'sum'],
                                                    help='Function to apply to values across feature to get a single '
                                                                                                   'value per feature.')
    summaryparser.add_argument('--category_summary', default='mean', choices = ['mean', 'median', 'min', 'max', 'sum'],
                                               help='Function to apply to feature values across category to get a single '
                                                                                                   'value per category.')
    summaryparser.add_argument('--chrom_len', type=argparse.FileType('r'), help="Tsv that provides with columns 'chrom' "
                                                                                "'length'. Not required, but considerably "
                                                                                "speeds up permutations.")

    topscoresparser.add_argument('--top_type',choices=['count','threshold','quantile'],
                                              help='How should the cutoff for top scores be determined?')
    topscoresparser.add_argument('--top', type=float, help='Cutoff for top scores in units determined by top_type.')
    g = topscoresparser.add_mutually_exclusive_group(required=False)
    g.add_argument("--ascending",dest="ascending", action="store_true", help="Sort ascending, (e.g., for p-values).")
    g.add_argument("--descending",dest="ascending", action="store_false", help="Sort descending, (e.g., for scores).")
    g.set_defaults(ascending=False)
    topscoresparser.add_argument('--max_dist', type=int, default=0,
                                                       help="Maximum distance (in bp) between feature and rod value "
                                                              "in order to consider the rod values for that feature.")


    mode_classes = {'Candidate':CandidateEnrichment,'Summary':SummaryEnrichment,
                    'TopScores':TopScoresEnrichment}

    def parse_cols(cols):
        cols1 = []
        for col in cols:
            try:
                cols1.append(int(col))
            except ValueError:
                cols1.append(col)
        return cols1



    args, unknown = parser.parse_known_args()

    help_str = "\n"+parser.format_help()
    if args.help:
        logger.setLevel(logging.INFO)



    if  args.run_type == 'Permute':
        permute_args, unknown = permuteparser.parse_known_args(unknown)
        if permute_args.mode is not None:
            modeparser = modeparsers.choices[permute_args.mode]
        else:
            permuteparser.error("argument --mode/-M is required")
        mode_args, unknown = modeparser.parse_known_args(unknown)
        help_str += "\n\n-----------Help for Permute step ------------\n\n" \
                                                           + permuteparser.format_help()
        help_str += "\n\n-----------Help for mode {} ------------\n\n".format(permute_args.mode) \
                                                           + modeparser.format_help()
    if args.run_type == 'Reduce':
        reduce_args, unknown = reduceparser.parse_known_args(unknown)
        help_str += "\n\n-----------Help for Reduce step ------------\n\n" \
                                                           + reduceparser.format_help()

    if args.help:
       logging.info(help_str)
       sys.exit(0)

    if unknown:
        logging.error("Unknown command line arguments: {}".format(unknown))
        logging.info(help_str)
        sys.exit(1)


    if args.feature_to_category is None:
        parser.error("argument --feature_to_category is required")

    logger.setLevel(getattr(logging,args.logging_level))


    logging.info("Loading feature to category mapping from {}.".format(args.feature_to_category.name))
    feature_to_category_cols = parse_cols(args.feature_to_category_cols)
    feature_to_category = pd.read_csv(args.feature_to_category,usecols=feature_to_category_cols,
                                                         sep=get_sep(args.feature_to_category.name))
    feature_to_category = feature_to_category[feature_to_category_cols]

    #test whether name is writeable:
    #we do not test all derived filenames, but this should be ok for most cases)
    #try:
    #    open(args.name,'w')
    #except Exception, e:
    #    logging.error("Cannot write to {}".format(args.name))
    #    raise
    dir0 = os.path.dirname(os.path.realpath(args.name))
    assert os.access(dir0, os.W_OK), \
                                "{} cannot be accessed for writing".format(dir0)

    if args.run_type == 'Permute':

        Enrichment = mode_classes[permute_args.mode]
        mode_args = {arg:getattr(mode_args,arg) for arg in vars(mode_args)}
        mode_args['ncpus'] = permute_args.ncpus

        if permute_args.mode in ['Summary', 'TopScores']:

            if mode_args['rod'] is None:
                modeparsers.choices[permute_args.mode].error("argument --rod is required")
            if mode_args['features'] is None:
                modeparsers.choices[permute_args.mode].error("argument --features is required")

            if permute_args.mode == 'Summary':
                if mode_args['chrom_len'] is not None:
                    chrom_len = pd.read_csv(mode_args['chrom_len'],sep=get_sep(mode_args['chrom_len'].name), squeeze=True,index_col=0)
                    mode_args['chrom_len'] = chrom_len


            logging.info("Loading features from {}.".format(mode_args['features'].name))
            feature_cols = parse_cols(mode_args.pop('feature_cols'))
            
            features_fh = mode_args.pop('features')
            feature_df = pd.read_csv(features_fh, index_col=[0,1],
                                        usecols=feature_cols, sep=get_sep(features_fh))
            feature_df.index.set_names(['chrom','start'], inplace=True)
            feature_df = feature_df[feature_cols[2:]]
            feature_df.rename(columns={feature_cols[2]:'end'}, inplace=True)

            if feature_df.columns[1] != feature_to_category.columns[0]:
                logging.warning("Feature name in features and feature_to_category do not match. "
                                 "They are {} and {}. Assuming that the feature identifiers in these "
                                 "columns are the same and using the first name. "
                                 "Make sure that these columns contain the same identifiers.".format(feature_df.columns[1],
                                                                                                    feature_to_category.columns[0]))
                feature_to_category.rename(columns={feature_to_category.columns[0]:feature_df.columns[1]}, inplace=True)


            logging.info("Loading rod from {}.".format(mode_args['rod'].name))
            rod_cols = parse_cols(mode_args.pop('rod_cols'))
            rod_fh = mode_args.pop('rod')
            rod_s = pd.read_csv(rod_fh, index_col=[0,1],
                                        usecols=rod_cols, sep=get_sep(rod_fh), squeeze=True)#,
                                        #header=False if type(rod_cols[0])==int else True)
            rod_s.index.set_names(['chrom','pos'], inplace=True)
            #remove inf values from input
            rod_s.replace([np.inf, -np.inf], np.nan, inplace=True)
            enrich = Enrichment(value_s=rod_s, feature_df=feature_df,
                                feature_to_category=feature_to_category, feature_name=feature_df.columns[1],
                                                               category_name=feature_to_category.columns[1], **mode_args)
        elif permute_args.mode == 'Candidate':
            if mode_args['candidates'] is None:
                modeparsers.choices[permute_args.mode].error("argument --candidates is required")
            candidate_features = []
            for candidate in mode_args.pop('candidates'):
                candidate_features.append(candidate.strip())
            enrich = Enrichment(candidate_features,feature_to_category=feature_to_category,
                                                               feature_name=feature_to_category.columns[0],
                                                               category_name=feature_to_category.columns[1], **mode_args)
        if not permute_args.noinfo:
            save_info(enrich, args.name)
            #enrich.create_info()
            #enrich.summary_per_feature.to_csv(args.name+'.summary_per_feature.tsv',sep='\t')
            #if permute_args.mode == 'TopScores':
            #    enrich.top_peaks.to_csv(args.name+'.top_peaks.tsv',sep='\t')
            #    enrich.peaks_per_gene.to_csv(args.name+'.peaks_per_gene.tsv',sep='\t')
        if permute_args.n_permut>0:
            start = time.time()
            enrich.permute(permute_args.n_permut)
            end = time.time()
            delta = end - start
            logging.info("{} permutations took took {} seconds = {} minutes = {} hours.".format(permute_args.n_permut,
                                                                                                delta,delta/60.,delta/3600.))

            if permute_args.noinfo:
                enrich.rank_table.to_csv(args.name+'.ranktable.tsv', sep='\t', header=True)
            else:
                enrich.get_pvals().to_csv(args.name+'.pvals.tsv', sep='\t', header=True)

    elif args.run_type == 'Reduce':
        if reduce_args.permuts is None:
            reduceparser.error("argument --permuts is required")

        permut_fhs = open_reduce_fns(reduce_args.permuts)
        tot_rank = reduce_fhs(permut_fhs)

        if reduce_args.category_to_description is not None:
            ctd_cols = parse_cols(reduce_args.category_to_description_cols)
            cat_to_desc = pd.read_csv(reduce_args.category_to_description,
                                                 usecols=ctd_cols, sep=get_sep(reduce_args.category_to_description))
        else:
            cat_to_desc = None

        p_val_df = rank_to_pval(tot_rank, pval_threshold=1,
                            category_to_description=cat_to_desc)
        p_val_df.index.name = 'category'
        #make this as a function to also use it in a all-in-one run

        try:
            logging.info("Loading {}".format(args.name+".peaks_per_feature.tsv"))
            cand_genes = np.unique(pd.read_csv(args.name+".peaks_per_feature.tsv",usecols=[0],sep='\t').values)
            #CONTINUE HERE.....
            gene_per_go_s = feature_to_category.set_index(feature_to_category_cols[0]).ix[cand_genes].groupby(feature_to_category_cols[1]).apply(lambda x: list(x.index))
            gene_per_go_s.name = feature_to_category_cols[0]
            p_val_df = pd.concat([p_val_df, gene_per_go_s], axis=1)
            def check_len(el):
                try:
                    return(len(el))
                except TypeError:
                    return 0
            #if not (p_val_df[feature_to_category_cols[0]].apply(check_len) == p_val_df["n_"+feature_to_category_cols[0]]).all():
            #    assert_df = p_val_df[["n_genes","genes"]]
            #    assert_df["len_genes"] = assert_df["genes"].apply(check_len)
            #    assert_df = assert_df[["n_genes","len_genes","genes"]]
            #    printv("Genes per category from",args.peaks_per_gene_fn,
            #            "inconsistent with n_genes reported in",permut_fhs[0].name,":",
            #            assert_df[assert_df["n_genes"]!=assert_df["len_genes"]])
        except IOError, e:
            logging.warning("Not adding feature names to pvalue file, no peaks_per_feature file found.")
            #raise e
            #logging.info(str(e))

        #try:
        #    p_val_df = p_val_df['n_genes']
        p_val_df.sort('p_value',inplace=True)
        p_val_df.to_csv(args.name+'_{}.pvals.tsv'.format(p_val_df['out_of'].max()), sep='\t')
        if reduce_args.remove_input:
            for fh in permut_fhs:
                os.remove(fh.name)


