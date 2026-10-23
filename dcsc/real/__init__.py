"""
Real-data validation of DCSC (see REAL_RESULTS.md):

    census_fetch  CELLxGENE Census subset (8 blood datasets, healthy vs COVID-19 / SLE)
    kang          IFN-β calibration contrast (Kang 2018, GSE96583)
    census        identity tree, IFN program, batch leakage, small panels (leave-one-dataset-out)
    norman        residual pursuit on real CRISPRa double perturbations (Norman 2019)
    sciplex       cross-cell-line drug Δ retrieval and dose read-out (sci-Plex3)
    twostage_eval two-stage retrieval on Norman, sci-Plex3 and Census cells
    census_1m     hybrid retrieval on a real 1,000,000-cell Census library
    donor_search  donor-level disease-state retrieval across labs (1M library + 3 external labs)
    donor_novelty test of the "no similar diseased donor" call (leave-one-collection-out)

Data are downloaded into data_cache/ (git-ignored).
"""
