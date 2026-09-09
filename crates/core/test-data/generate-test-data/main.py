import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csc_matrix, csr_matrix

# First 100 genes from ../../src/gene_list/chemistry/xenium_v1_human.rs
GENES = [
    ("ENSG00000116678", "LEPR"),
    ("ENSG00000258377", "AL139099.1"),
    ("ENSG00000267057", "LINC01905"),
    ("ENSG00000265413", "AP001094.2"),
    ("ENSG00000236572", "AC006369.2"),
    ("ENSG00000240889", "NDUFB2-AS1"),
    ("ENSG00000120802", "TMPO"),
    ("ENSG00000285748", "AC090337.2"),
    ("ENSG00000106459", "NRF1"),
    ("ENSG00000010017", "RANBP9"),
    ("ENSG00000277945", "AC107308.1"),
    ("ENSG00000189366", "ALG1L"),
    ("ENSG00000105088", "OLFM2"),
    ("ENSG00000235615", "AJ239322.1"),
    ("ENSG00000198908", "BHLHB9"),
    ("ENSG00000223774", "AL513217.1"),
    ("ENSG00000284988", "LINC02203"),
    ("ENSG00000125788", "DEFB126"),
    ("ENSG00000172757", "CFL1"),
    ("ENSG00000120054", "CPN1"),
    ("ENSG00000286803", "AL109824.1"),
    ("ENSG00000238074", "TSPY9P"),
    ("ENSG00000184164", "CRELD2"),
    ("ENSG00000285738", "AC098969.1"),
    ("ENSG00000286329", "AC108690.1"),
    ("ENSG00000257105", "AC006581.1"),
    ("ENSG00000263435", "AC024610.1"),
    ("ENSG00000230651", "RGPD4-AS1"),
    ("ENSG00000232520", "GCSIR"),
    ("ENSG00000110851", "PRDM4"),
    ("ENSG00000144488", "ESPNL"),
    ("ENSG00000287150", "AL139340.1"),
    ("ENSG00000011677", "GABRA3"),
    ("ENSG00000177725", "AC105206.1"),
    ("ENSG00000048471", "SNX29"),
    ("ENSG00000278467", "AC138393.3"),
    ("ENSG00000236372", "AL353052.2"),
    ("ENSG00000146809", "ASB15"),
    ("ENSG00000132676", "DAP3"),
    ("ENSG00000180475", "OR10Q1"),
    ("ENSG00000228340", "MIR646HG"),
    ("ENSG00000274031", "AC092140.2"),
    ("ENSG00000249001", "AC093895.1"),
    ("ENSG00000285895", "AP003557.2"),
    ("ENSG00000185467", "KPNA7"),
    ("ENSG00000284602", "AL031432.4"),
    ("ENSG00000228649", "SNHG26"),
    ("ENSG00000268112", "AC008761.1"),
    ("ENSG00000265786", "LINC01906"),
    ("ENSG00000110756", "HPS5"),
    ("ENSG00000235435", "AC096666.1"),
    ("ENSG00000259119", "AL132796.2"),
    ("ENSG00000259721", "AC090877.2"),
    ("ENSG00000177302", "TOP3A"),
    ("ENSG00000204033", "LRIT2"),
    ("ENSG00000150455", "TIRAP"),
    ("ENSG00000184814", "PRR23B"),
    ("ENSG00000117450", "PRDX1"),
    ("ENSG00000260563", "AC132872.1"),
    ("ENSG00000238007", "AC024619.1"),
    ("ENSG00000111331", "OAS3"),
    ("ENSG00000167081", "PBX3"),
    ("ENSG00000178358", "OR2D3"),
    ("ENSG00000132846", "ZBED3"),
    ("ENSG00000176593", "AC008969.1"),
    ("ENSG00000287068", "AP004606.1"),
    ("ENSG00000120215", "MLANA"),
    ("ENSG00000224614", "TNK2-AS1"),
    ("ENSG00000269437", "NXF2B"),
    ("ENSG00000286477", "AC095038.4"),
    ("ENSG00000273841", "TAF9"),
    ("ENSG00000168385", "SEPTIN2"),
    ("ENSG00000234597", "AC010096.1"),
    ("ENSG00000269973", "AC010969.2"),
    ("ENSG00000221818", "EBF2"),
    ("ENSG00000224712", "NPIPA3"),
    ("ENSG00000114107", "CEP70"),
    ("ENSG00000137338", "PGBD1"),
    ("ENSG00000287871", "AC073632.1"),
    ("ENSG00000198054", "DSCR8"),
    ("ENSG00000223850", "MYCNUT"),
    ("ENSG00000115216", "NRBP1"),
    ("ENSG00000135373", "EHF"),
    ("ENSG00000198570", "RD3"),
    ("ENSG00000256196", "AP003721.1"),
    ("ENSG00000197408", "CYP2B6"),
    ("ENSG00000269113", "TRABD2B"),
    ("ENSG00000249898", "MCPH1-AS1"),
    ("ENSG00000121316", "PLBD1"),
    ("ENSG00000145632", "PLK2"),
    ("ENSG00000283095", "FP565171.1"),
    ("ENSG00000258753", "AL359792.1"),
    ("ENSG00000100373", "UPK3A"),
    ("ENSG00000089123", "TASP1"),
    ("ENSG00000105617", "LENG1"),
    ("ENSG00000253821", "AC090796.1"),
    ("ENSG00000147180", "ZNF711"),
    ("ENSG00000227155", "AL161725.1"),
    ("ENSG00000162302", "RPS6KA4"),
    ("ENSG00000230732", "AC016949.1"),
]


def main():
    # Adapted from https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html
    ensembl_ids, gene_names = zip(*GENES)
    n_cells = 10

    rng = np.random.default_rng()
    counts = rng.integers(0, 10, size=(n_cells, len(GENES)))
    counts[0] = 10

    csr_adata = ad.AnnData(csr_matrix(counts, dtype=np.float32))
    csc_adata = ad.AnnData(csc_matrix(counts, dtype=np.float32))
    dense_adata = ad.AnnData(counts)

    for name, adata in [
        ("csr_adata", csr_adata),
        ("csc_adata", csc_adata),
        ("dense_adata", dense_adata),
    ]:
        adata.layers["counts"] = adata.X.copy()
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.obs["barcode"] = adata.obs_names
        adata.obs["annotation"] = [f"group{i % 2}" for i in range(len(adata.obs.index))]

        unannotated_cell = list(adata.obs["annotation"])[:-1] + [None]
        adata.obs["annotation_missing"] = unannotated_cell

        adata.var_names = list(gene_names)
        adata.var["ensembl_id"] = list(ensembl_ids)
        adata.var["feature_types"] = ["Gene Expression"] * len(adata.var)

        adata.write_h5ad(f"../{name}.h5ad")

    filename = "../1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5"
    tenx_adata = sc.read_10x_h5(filename)

    tenx_adata.obs["barcode"] = tenx_adata.obs.index
    tenx_adata.obs["annotation"] = [
        f"group_{i % 2}" for i in range(len(tenx_adata.obs.index))
    ]
    tenx_adata.var["gene_name"] = tenx_adata.var_names

    sc.write(f"{filename}ad", tenx_adata)


if __name__ == "__main__":
    main()
