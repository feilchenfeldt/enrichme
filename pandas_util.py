import numpy as np
import pandas as pd
"""
Specifications (So far, only implemented for the single index part below):


feature_df ... Data Frame of intervals along the genome,
                equivalent of a bed file, but 1-indexed
index: 
    (chrom, start)
columns 
    required: 
        end
    optional:
        name
        ...
"""

def check_rod(rod):
    """
    Check if input follows the ROD specifications.
    Throw assertion error otherwise.
    Specifications:
    rod ... Series or DataFrame of reference ordered data
    index:
        (chrom, pos)
    """
    assert rod.index.names == ('chrom','pos'), "Index names must be ('chrom','pos')."
    assert rod.ix and rod.iloc and rod.loc, ("ROD lacks pandas indexing functionality. "
                                                            "Is it DataFrame or Series?")

def check_feature_df(feature_df):
    """
    Check if input follows the ROD specifications.
    Throw assertion error otherwise.
    Specifications:
    index:
        (chrom, start)
    columns
        required:
            end
        optional:
            feature
    """
    assert feature_df.index.names == ('chrom','start'), "Index names must be ('chrom','start')."
    assert feature_df.ix and feature_df.iloc and feature_df.loc, ("feature_df lacks pandas indexing "
                                                                  "functionality. Is it DataFrame or Series?")
    assert 'end' in feature_df.columns, "feature_df lacks required column 'end'"

def index_rolling(s,window,func,overlap=0,min_n_values=0,*args,**kwargs):
    """
    Apply function in rolling windows, where the window
    size is defined with respect to the index values.
    This means that different windows can comprise different
    numbers of elements.
    
    s ... pandas Series
    window ... window size in units of the index values
    func ... function to apply to the series values within each
                window
    overlap ... oberlap size of windows
    min_n_values ... minimum number of values (e.g. SNPs) per window
                    None will be returned for windows that have less values.

    args, kwarg ... additional arguments for func
    
    Example: index_rolling(pd.Series(range(20),index=np.random.rand(20))),0.1,max)
    """
    assert isinstance(s, pd.Series)
    #note that basis must be sorted in order for this to work properly
    windows_min = s.index.min()
    windows_max = s.index.max()
    window_starts = np.arange(windows_min, windows_max-window, window-overlap)
    window_starts = pd.Series(window_starts, index = window_starts+window/2)
    def applyToWindow(val):
        # using slice_indexer rather that what.loc [val:val+window] allows
        # window limits that are not specifically in the index
        try:
            indexer = s.index.slice_indexer(val,val+window,1)
        except IndexError:
            print val, val+window
            print s
            raise
        chunk = s.iloc[indexer]
        try:
            if len(chunk) < min_n_values:
                #print indexer, "too few snps"
                return None
        except TypeError:
            return None
        try:
            return func(chunk,*args,**kwargs)
        except ValueError, e:
            if "empty sequence" in str(e):
                #print indexer, chunk
                return None
            else:
                raise
    rolled = window_starts.apply(applyToWindow)
    return rolled


##How should the output of this look like? What index?

def data_per_feature(rod,feature_df, feature_name='feature', max_dist=0):
    """
    Get the entires in rod which lie within a feature
    (e.g. gene) in feature_df.
    Input:
    rod (reference ordered data)... pandas series or data frame with multiindex (chrom, pos)
                                    such as SNP genotypes
    feature_df (gene annotation data frame)... index must be (chrom,feature_name), must have columns 'start', 'end'
    
    max_dist not implemented yet, we need to take care of overlapping genes,
    currently we only take the last
    """
    rod = pd.DataFrame(rod)
    chrom_features = []
    chrpos_names = rod.index.names
    for chrom in rod.index.droplevel(1).unique():
        try:
            feature_chrom = feature_df.ix[chrom]
        except KeyError:
            continue
        rod_chrom = rod.ix[chrom]
        if not feature_chrom.index.is_monotonic:
            feature_chrom = feature_chrom.sort_index()
        if not rod_chrom.index.is_monotonic:
            rod_chrom = rod_chrom.sort_index()
        pos_rel_to_start = feature_chrom.index.searchsorted(rod_chrom.index)
        pos_rel_to_end = np.searchsorted(feature_chrom["end"].values,rod_chrom.index.values)
        in_feature = (pos_rel_to_start - pos_rel_to_end) == 1
        feature_id = feature_chrom.iloc[pos_rel_to_end[in_feature]][feature_name].values
        snp_df = rod_chrom[in_feature].copy()
        snp_df['chrom'] = chrom
        snp_df[feature_name] = feature_id
        chrom_features.append(snp_df)
    dpf = pd.concat(chrom_features)
    dpf.set_index(['chrom'],append=True,inplace=True)
    dpf = dpf.reorder_levels(['chrom','pos'])
    return dpf


def get_features_per_data(peak_s, feature_df, feature_name='feature', max_dist=0):
    """
    take the input data series and gets a similar series
    with one entry per pair data-point gene
    (i.e., there can be 0,1 or more entries per data point)
    
    """
    all_features = []
    if not feature_df.index.is_monotonic:
        feature_df = feature_df.sort_index()
    tot_hit_df = pd.DataFrame()
    for chrom in peak_s.index.droplevel(1).unique():
        loc_feature_df = feature_df.ix[chrom]
        #loc_feature_df = loc_feature_df.append(pd.DataFrame(np.nan,index=[np.inf],columns=loc_feature_df.columns))
        #print loc_feature_df.index-max_dist, peak_s.ix[chrom].index.values
        #try:
        pos_rel_to_start = np.searchsorted(loc_feature_df.index.values-max_dist,peak_s.ix[chrom].index.values)
        #except:
        #    print chrom, peak_s.ix[chrom]
        pos_rel_to_end = np.searchsorted(loc_feature_df["end"].values+max_dist,peak_s.ix[chrom].index.values)
        features_per_datapoint = (pos_rel_to_start - pos_rel_to_end)
        #why is this so slow print features_per_datapoint.shape
        data_idx = [i for i in range(len(features_per_datapoint)) for j in range(features_per_datapoint[i])]
        features = loc_feature_df[feature_name].iloc[np.hstack([range(a,b) for a,b in zip(pos_rel_to_end,pos_rel_to_start)])].values
        data_df = pd.DataFrame(peak_s.ix[chrom].iloc[data_idx])
        data_df[feature_name] = features
        data_df['chrom'] = chrom
        all_features.append(data_df)
    fpd=pd.concat(all_features)    
    #fpd.set_index(['chrom'],append=True,inplace=True)
    #fpd = dpf.reorder_levels(['chrom','pos'])
    g = fpd.reset_index().groupby(['chrom','pos'])
    def get_series_of_features(gdf):
        features = gdf[feature_name].unique()
        r = pd.Series({i:s for i,s in enumerate(features)})
        df = pd.DataFrame({feature_name:r,peak_s.name:gdf[peak_s.name].values[0]})
        return df
    d = g.apply(get_series_of_features)
    d.index.names = ['chrom','pos','number']
    return  d

def features_per_data_to_data_per_features(fpd, feature_name='features'):
    return fpd.reset_index().set_index(feature_name).sort_index()



def get_features(peak_s, feature_df, feature_name='feature', max_dist=0):
    """
    take the input series and gets.
    names of features nearby

    Input:
    peak_s ... pandas series with (chrom, pos) index and value of
                the statistic ('peak height'). Series should be named.
    feature_df ... data frame with feature info.
    """
    all_features = []
    if not feature_df.index.is_monotonic:
        feature_df = feature_df.sort_index()
    tot_hit_df = pd.DataFrame()
    for chrom in peak_s.index.droplevel(1).unique():
        loc_feature_df = feature_df.ix[chrom]
        #loc_feature_df = loc_feature_df.append(pd.DataFrame(np.nan,index=[np.inf],columns=loc_feature_df.columns))
        #print loc_feature_df.index-max_dist, peak_s.ix[chrom].index.values
        #try:
        pos_rel_to_start = np.searchsorted(loc_feature_df.index.values-max_dist,peak_s.ix[chrom].index.values)
        #except:
        #    print chrom, peak_s.ix[chrom]
        pos_rel_to_end = np.searchsorted(loc_feature_df["end"].values+max_dist,peak_s.ix[chrom].index.values)
        features = list(set(loc_feature_df[feature_name].iloc[np.hstack([range(a,b) for a,b in zip(pos_rel_to_end,pos_rel_to_start)])]))
        all_features += features
    return all_features


def apply_to_feature(feature_df,groupby_func_name=None,function=None):
    """
    Apply a function to the entries for each feature.
    feature_df ... dataframe with index (chrom, feature_name, pos)
                   (Such as the output of data_per_feature())
    groupby_func_name ... name of the function of the groupby object
                         to apply to the data
                          This is faster than applying a function object.
    function ... alternatively: function object to apply
    """
    groups = feature_df.groupby(lambda idx: idx[1])
    if groupby_func_name is not None:
        return getattr(groups,groupby_func_name)()
    elif function is not None:
        return groups.apply(function)
    else:
        raise ValueError("Either groupby_func_name or function have to be given.")


#--------------WORK WITH SINGLE NUMERIC INDEX------------------------


def rod_to_1d(rod, chrom_len_s, drop=True):
    """
    Converts a (chrom, pos) multiindex to a 
    single numeric index that runs through chromosomes.
    Note that the (chrom, pos) index is sorted lexographically,
    i.e., if chrom names are strings, the entries in the resulting
    index are Chr1, Chr10, Chr2, Chr20, Chr3,...,ChrX, ChrY.
    Example: (Chr2, 1) is converted to len(Chr1)+len(Chr10)+1.
    The inverse operation is given by rod_to_chrompos.
    
    Input:
    rod ... series or dataframe with reference ordered data 
            with multiindex (chrom,pos)
    chrom_len_s ... series with chromosome names as keys
                    and chromosome length as values
    drop ... If False, keep chrom, pos as columns
    """
    if not chrom_len_s.index.is_monotonic:
        chrom_len_s = chrom_len_s.sort_index()
    rod = pd.DataFrame(rod).copy()
    if not rod.index.is_monotonic:
        rod.sortlevel(inplace=True)
    columns = columns = [c for c in rod.columns if c not in ['index']]
    rod.reset_index(inplace=True)
    #return rod.groupby('chrom')
    try:
        index = rod.groupby('chrom').apply(lambda df: df['pos']+chrom_len_s.loc[:df['chrom'].iloc[0]].iloc[:-1].sum()).values
    except KeyError, e:
        print chrom_len_s
        raise e
    rod['index'] = index
    rod.set_index('index', inplace=True, drop=True)
    if not drop:
        columns = list(columns) + ['chrom', 'pos']
    rod.index = rod.index.values.astype(int)
    if not rod.index.is_monotonic:
        rod.sort_index(inplace=True)
    return rod[columns] if len(columns)>1 else rod[columns[0]]

def rod_to_chrompos(rod_1d, chrom_len_s, drop=True):
    """
    Reverts the action of rod_to_1d.
    Converts a single numeric index that runs through 
    chromosomes to a (chrom, pos) multiindex.
    Note that the single index is expected to correspond
    to (chrom, pos) sorted lexographically,
    i.e., if chrom names are strings, the entries should be in
    in the order Chr1, Chr10, Chr2, Chr20, Chr3,...,ChrX, ChrY.
    Example:  len(Chr1)+len(Chr10)+1 is converted to (Chr2, 1).
    
    Input:
    rod ... series or dataframe of reference ordered data 
            with single index running through chromosomes.
            (as produced by rod_to_1d())
    chrom_len_s ... series with chromosome names as keys
                    and chromosome length as values
    drop ... If False, keep numeric index as columns
    """
    if not chrom_len_s.index.is_monotonic:
        chrom_len_s = chrom_len_s.sort_index()
    rod = pd.DataFrame(rod_1d).copy()
    if not rod.index.is_monotonic:
        rod.sort_index(inplace=True)
    columns = [c for c in rod.columns if c not in ['chrom','pos']]
    cs = chrom_len_s.cumsum()
    rod['chrom'] = np.nan
    for chrom, (start, end) in  zip(cs.index,zip([0] + list(cs.values[:-1]),cs.values)):
        end = min(end,rod.index[-1])
        rod.loc[slice(start,end),'chrom'] = chrom
        #print chrom, rod.loc[slice(start,end),columns[0]]
        rod.loc[slice(star/t,end),'pos'] = rod.ix[slice(start,end)].index - start
    if drop:
        rod.set_index(['chrom','pos'], inplace=True, drop=True)
    else:
        rod = rod.reset_index().set_index(['chrom','pos'], drop=True)
        columns = list(columns) + ['index']
    #if not rod.index.is_monotonic:
    #    rod.sort_index(inplace=True)

    return rod[columns] if len(columns)>1 else rod[columns[0]]

def feature_df_to_1d(feature_df, chrom_len_s):
    """
    Converts mulitindex feature_df (chrom,start)
    to single numeric index running trough all 
    chromosomes. The column 'end' is also converted.
    See rod_to_1d for details.
    """
    feature_df = feature_df.copy()
    feature_df.index.names = (feature_df.index.names[0], 'pos') 
    feature_df_1d = rod_to_1d(feature_df, chrom_len_s,drop=False)
    #print feature_df_1d
    end_df = feature_df_1d.reset_index().set_index(['chrom','end'])
    #print end_df
    end_df.drop('pos',axis=1,inplace=True)
    end_df.rename(columns={'index':'start'}, inplace=True)
    end_df.index.names = (end_df.index.names[0], 'pos') 
    end_1d = rod_to_1d(end_df, chrom_len_s)
    end_1d.index.name = 'end'
    end_1d = end_1d.reset_index().set_index('start')
    end_1d.index.name = 'index'
    if not end_1d.index.is_monotonic:
        end_1d.sort_index(inplace=True)
    return end_1d

def feature_df_to_chrompos(feature_df_1d, chrom_len_s):
    """
    Converts feature_df with single numberic index
    running through all chromosomes to 
    multiindex (chrom,start)
    The column 'end' is also converted.
    This is the inverse function of 
    feature_df_to_1d().
    See rod_to_chrompos for details.
    """
    feature_df_1d = feature_df_1d.copy()
    feature_df = rod_to_chrompos(feature_df_1d, chrom_len_s)
    feature_df.index.names = (feature_df.index.names[0], 'start') 
    end_df_1d = feature_df.reset_index().set_index('end')
    end_df = rod_to_chrompos(end_df_1d, chrom_len_s)
    end_df.index.names = (end_df.index.names[0],'end')
    end_df = end_df.reset_index().set_index(['chrom','start'])
    if not end_df.index.is_monotonic:
        end_df.sortlevel(inplace=True)
    return end_df

def data_per_feature_FI(rod, feature_df, feature_name = 'feature'):
    """
    Get the entires in rod which lie within a feature
    (e.g. gene) in feature_df. FI stands for flattened index.
    Input:
    rod (reference ordered data)... pandas series or data frame with numeric index
                                    such as SNP genotypes
    feature_df (gene annotation data frame)... index must be (chrom,feature_name), must have columns 'start', 'end'
    """
    rod = pd.DataFrame(rod)
    pos_rel_to_start = feature_df.index.searchsorted(rod.index)
    #feature_df["end"] is not necessarily sorted, but wouldn't sorting it here
    #lead to problems  as well?
    pos_rel_to_end = np.searchsorted(feature_df["end"].values,rod.index.values)
    in_feature = (pos_rel_to_start - pos_rel_to_end) == 1
    feature_id = feature_df.iloc[pos_rel_to_end[in_feature]][feature_name].values
    rod = rod[in_feature].copy()
    rod[feature_name] = feature_id
