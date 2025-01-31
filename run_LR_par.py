import numpy as np
import os
import scanpy as sc
import pandas as pd
import scipy.spatial as spatial
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=int)
parser.add_argument("-l", "--LRs_file", type=str)
parser.add_argument("-p", "--pucks_file", type=str)
parser.add_argument("-r", "--results_dir", type=str)
parser.add_argument("-n", "--n_pairs", type=int)
parser.add_argument("-m", "--m_pucks", type=int)

args = parser.parse_args()

cores = args.cores
LRs_file = args.LRs_file
pucks_file = args.pucks_file
results_dir = args.results_dir
n_pairs = args.n_pairs
m_pucks = args.m_pucks


def select_LR_pairs(LRs_file, n_pairs):
    LRs = pd.read_csv(LRs_file, index_col=0)

    LRs["pair"] = LRs["ligand.complex"] + "--" + LRs["receptor.complex"]

    pairs = LRs["pair"].unique()

    if n_pairs is None:
        return pairs
    else:
        return pairs[:n_pairs]


def select_pucks_df(pucks_file, m_pucks):
    pucks_df = pd.read_csv(pucks_file)

    if m_pucks is None:
        return pucks_df
    else:
        return pucks_df[:m_pucks]


def LR_calcs(pair, pucks_df, results_dir):
    print(pair, "started")

    if not os.path.exists(results_dir):
        try:
            os.makedirs(results_dir)
            print(f"{results_dir} created")
        except Exception as e:
            print(f"A {e} occurred, {results_dir} have already been created")
    else:
        print(f"{results_dir} already exists")

    ligand_centered = True
    N = 100
    results = []
    discarded_LR = []
    # for puck in test_pucks:
    for pp in range(len(pucks_df)):
        puck = pucks_df.iloc[pp].puck

        adata_path = pucks_df.iloc[pp].adata_path

        patient = pucks_df.iloc[pp].patient

        print(puck, adata_path, patient)

        ligand = pair.split("--")[0].split("_")
        receptor = pair.split("--")[1].split("_")

        adata = sc.read(adata_path)

        # print(adata.shape)
        bead_thresh = np.quantile(adata.obs.n_genes_by_counts, 0.25)
        bead_thresh = np.max([bead_thresh, 150])
        sc.pp.filter_cells(adata, min_genes=bead_thresh)
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)

        # fig, ax = plt.subplots(nrow, ncol, figsize=(ncol*wd,nrow*wd), gridspec_kw={'wspace':0.2})
        # axr=ax.ravel()

        puck_locs = adata.obsm["spatial"]

        # ligand_counts=adata[(adata[:,ligand].X>0),ligand].X.todense()
        # receptor_counts=adata[(adata[:,receptor].X>0),receptor].X.todense()

        if len(ligand) > 1:
            try:
                shares = list(set(adata.var.index) & set(ligand))
                ad = adata[:, shares].copy()
                # sc.pp.filter_cells(ad, min_genes=1)
                sc.pp.filter_cells(ad, min_counts=1)
                qun = np.quantile(ad.obs.n_counts, 0.75)
                sc.pp.filter_cells(ad, min_counts=qun)
                ligand_locs = ad.obsm["spatial"]
            except:
                ligand_locs = np.array([])

        else:
            try:
                gene_counts = (
                    adata[adata[:, ligand].X > 0][:, ligand].X.todense().tolist()
                )
                thresh = np.quantile(gene_counts, 0.75)
                ligand_locs = adata[(adata[:, ligand].X >= thresh)].obsm["spatial"]
            except:
                ligand_locs = np.array([])

        if len(receptor) > 1:
            try:
                shares = list(set(adata.var.index) & set(receptor))
                ad = adata[:, shares].copy()
                # sc.pp.filter_cells(ad, min_genes=1)
                sc.pp.filter_cells(ad, min_counts=1)
                qun = np.quantile(ad.obs.n_counts, 0.75)
                sc.pp.filter_cells(ad, min_counts=qun)
                receptor_locs = ad.obsm["spatial"]
            except:
                receptor_locs = np.array([])

        else:
            try:
                gene_counts = (
                    adata[adata[:, receptor].X > 0][:, receptor].X.todense().tolist()
                )
                thresh = np.quantile(gene_counts, 0.75)
                receptor_locs = adata[(adata[:, receptor].X >= thresh)].obsm["spatial"]
            except:
                receptor_locs = np.array([])

        # tumor_tcr_locs=ad_tcr[ad_tcr.obs.tumor_beads==1].obsm['spatial']
        # print(pair,puck,receptor_locs.shape[0],ligand_locs.shape[0])

        if receptor_locs.shape[0] > 20 and ligand_locs.shape[0] > 20:
            # to speed up we can switch these
            if ligand_locs.shape[0] < receptor_locs.shape[0]:
                ligand_centered = True
            else:
                ligand_centered = False

            for j, d in enumerate([15, 30, 100, 300]):
                if ligand_centered:
                    # print('ligand_centered')
                    point_tree = spatial.cKDTree(receptor_locs)
                    t_N = point_tree.query_ball_point(ligand_locs, d)
                    true = sum([len(n) for n in t_N])
                    counts = []
                    for i in range(N):
                        fake_ligands = puck_locs[
                            np.random.choice(
                                puck_locs.shape[0],
                                size=ligand_locs.shape[0],
                                replace=False,
                            ),
                            :,
                        ]
                        t_N = point_tree.query_ball_point(fake_ligands, d)
                        counts.append(sum([len(n) for n in t_N]))
                else:
                    # print('receptor_centered')
                    point_tree = spatial.cKDTree(ligand_locs)
                    t_N = point_tree.query_ball_point(receptor_locs, d)
                    true = sum([len(n) for n in t_N])
                    counts = []
                    for i in range(N):
                        fake_receptors = puck_locs[
                            np.random.choice(
                                puck_locs.shape[0],
                                size=receptor_locs.shape[0],
                                replace=False,
                            ),
                            :,
                        ]
                        t_N = point_tree.query_ball_point(fake_receptors, d)
                        counts.append(sum([len(n) for n in t_N]))

                counts = np.array(counts)
                # min_c=counts.min();max_c=counts.max()
                # y, x, _ = axr[f].hist(counts, range=(min_c,max_c), bins=int(max_c-min_c)+1)

                # print(sum(counts),true)
                pval = np.round((sum(counts > true) + 1) / (N + 1), 7)

                # plt.title('Sum of number of ERV beads at '+str(d)+' pixels radius of each TCR')
                # plt.title('Sum of all ERV exressing beads \n
                #          in 'str(d)+' pixels radius of each TCR vs random shuffles of N locations from')
                # axr[f].set_title(pair + '  ' + '\n'+str(d)+' pixels radius'+ '\n'+'p = '+str(pval) )
                # axr[f].axvline(true, 0, ymax=y.max(), color='r')

                results.append(
                    [
                        pair,
                        puck,
                        patient,
                        pval,
                        true,
                        np.median(counts),
                        np.mean(counts),
                        np.std(counts),
                        d,
                        receptor_locs.shape[0],
                        ligand_locs.shape[0],
                    ]
                )
        else:
            discarded_LR.append(
                [pair, puck, receptor_locs.shape[0], ligand_locs.shape[0]]
            )

    try:
        res = pd.DataFrame(results)
        res.columns = [
            "LR",
            "puck",
            "patient",
            "pval",
            "hits",
            "median",
            "mean",
            "std",
            "distance",
            "receptor_beads",
            "ligand_beads",
        ]
        # print(res)
        res.to_csv(f"{results_dir}/result_{pair}.csv")
    except:
        print(pair, "finished")
        # print(f'pair {pair} has no sample')
    try:
        dis = pd.DataFrame(discarded_LR)
        dis.columns = ["LR", "puck", "receptor_beads", "ligand_beads"]
        dis.to_csv(f"{results_dir}/discard_{pair}.csv")
    except:
        print(pair, "finished")
        # print(f'pair {pair} has all samples')


pairs = select_LR_pairs(LRs_file, n_pairs)
pucks_df = select_pucks_df(pucks_file, m_pucks)
args = [(pair, pucks_df, results_dir) for pair in pairs]

# print(args)

pool = Pool(cores)
results = pool.starmap(LR_calcs, args)
pool.close()
