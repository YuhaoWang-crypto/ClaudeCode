-- Cell state library v0.2: long-format schema.
-- One row per (reference profile, gene). A profile is one cell type x one condition
-- measured against a named control. Adding a cell type, stimulus, dose or time point
-- is a new profile, not a new column.
PRAGMA foreign_keys = ON;

CREATE TABLE sources(
  source_id      TEXT PRIMARY KEY,
  url            TEXT NOT NULL,
  version        TEXT,              -- e.g. git commit
  sha256         TEXT,
  evidence_type  TEXT NOT NULL,     -- published_aggregate_response_matrix | own_measurement | ...
  access_status  TEXT NOT NULL,
  license_note   TEXT
);

-- Cell types form a tree (PBMC > L1 > L2). Author labels map in via cell_type_aliases.
CREATE TABLE cell_types(
  cell_type_id   TEXT PRIMARY KEY,
  display_name   TEXT NOT NULL,
  lineage_level  TEXT NOT NULL,     -- bulk | L1 | L2
  parent_id      TEXT REFERENCES cell_types(cell_type_id),
  note           TEXT
);

CREATE TABLE cell_type_aliases(
  source_id      TEXT NOT NULL REFERENCES sources(source_id),
  alias_level    TEXT NOT NULL,     -- label level used by the source (L1 / L2 / PBMC)
  alias_label    TEXT NOT NULL,
  cell_type_id   TEXT NOT NULL REFERENCES cell_types(cell_type_id),
  PRIMARY KEY(source_id, alias_level, alias_label)
);

-- Stimulus x dose x time. Cell preparation lives on the profile, not here.
CREATE TABLE conditions(
  condition_id   TEXT PRIMARY KEY,
  stimulus       TEXT NOT NULL,
  dose           REAL,
  dose_unit      TEXT,
  duration_h     REAL,
  source_url     TEXT,
  status         TEXT NOT NULL
);

CREATE TABLE reference_profiles(
  profile_id           TEXT PRIMARY KEY,
  source_id            TEXT NOT NULL REFERENCES sources(source_id),
  cell_type_id         TEXT NOT NULL REFERENCES cell_types(cell_type_id),
  condition_id         TEXT NOT NULL REFERENCES conditions(condition_id),
  control_condition_id TEXT NOT NULL REFERENCES conditions(condition_id),
  culture_context      TEXT NOT NULL,  -- how the cells were prepared when stimulated
  n_donors             INTEGER,
  summary_stat         TEXT NOT NULL,  -- what one value means (e.g. median log2FC across donors)
  gene_selection       TEXT NOT NULL,  -- why these genes are in the profile
  status               TEXT NOT NULL,  -- active | active_partial | staged_parent_level | ...
  UNIQUE(source_id, cell_type_id, condition_id, control_condition_id)
);

CREATE TABLE reference_values(
  profile_id  TEXT NOT NULL REFERENCES reference_profiles(profile_id),
  gene        TEXT NOT NULL,
  log2fc      REAL NOT NULL,
  se          REAL,                 -- NULL: not available in the aggregate source
  PRIMARY KEY(profile_id, gene)
);
CREATE INDEX idx_reference_values_gene ON reference_values(gene);

-- Per cell type: can the stimuli be told apart on this cell type's genes?
CREATE TABLE identifiability(
  cell_type_id      TEXT PRIMARY KEY REFERENCES cell_types(cell_type_id),
  assumed_sigma     REAL NOT NULL,  -- assumed per-gene noise of a target log2FC (assumption, not measured)
  n_genes           INTEGER NOT NULL,
  rms_ifna          REAL NOT NULL,
  rms_ifng          REAL NOT NULL,
  cosine            REAL NOT NULL,
  condition_number  REAL NOT NULL,
  se_ifna           REAL NOT NULL,  -- SE of the coefficient, in units of one full reference response
  se_ifng           REAL NOT NULL,
  z_ifna            REAL NOT NULL,  -- 1/SE: z-score a full-strength response would reach
  z_ifng            REAL NOT NULL,
  status            TEXT NOT NULL,  -- both | ifna_only | ifng_only | none
  author_ifng_zero_rule INTEGER NOT NULL, -- 1 if the author R code would force IFNg_score to 0
  notes             TEXT
);

-- Gene-pair code panels (one per cell type). Each pair is a ternary bit: sign(b_i - b_j) with dead zone tau.
CREATE TABLE pair_panels(
  panel_id       TEXT NOT NULL,
  cell_type_id   TEXT NOT NULL REFERENCES cell_types(cell_type_id),
  rank           INTEGER NOT NULL,
  gene_i         TEXT NOT NULL,
  gene_j         TEXT NOT NULL,
  d_ifna         REAL NOT NULL,     -- reference IFNa log2FC(i) - log2FC(j)
  d_ifng         REAL NOT NULL,
  tau            REAL NOT NULL,
  sigma          REAL NOT NULL,
  PRIMARY KEY(panel_id, rank)
);

-- Own measurements (still empty in v0.2) and planned experiments, carried over from v0.1.
CREATE TABLE samples(
  sample_id         TEXT PRIMARY KEY,
  donor_id          TEXT,
  cell_type_id      TEXT REFERENCES cell_types(cell_type_id),
  condition_id      TEXT REFERENCES conditions(condition_id),
  paired_control_id TEXT REFERENCES samples(sample_id),
  culture_context   TEXT,
  protocol_version  TEXT,
  assay             TEXT,
  batch_id          TEXT,
  replicate_kind    TEXT,
  time_h            REAL,
  viability_pct     REAL,
  data_uri          TEXT,
  status            TEXT
);

CREATE TABLE observations(
  sample_id   TEXT REFERENCES samples(sample_id),
  feature     TEXT,
  measurement REAL,
  unit        TEXT,
  modality    TEXT,
  PRIMARY KEY(sample_id, feature, modality)
);

CREATE TABLE planned_experiments(
  plan_id         TEXT PRIMARY KEY,
  donor_slot      TEXT,
  cell_type_id    TEXT REFERENCES cell_types(cell_type_id),
  condition_id    TEXT REFERENCES conditions(condition_id),
  control_plan_id TEXT REFERENCES planned_experiments(plan_id),
  status          TEXT
);

-- Wide view for anyone who wants the v0.1 shape back (one row per cell type x gene).
CREATE VIEW v_reference_wide AS
SELECT p.cell_type_id, v.gene,
       MAX(CASE WHEN p.condition_id = 'IFNA' THEN v.log2fc END) AS ifna_log2fc,
       MAX(CASE WHEN p.condition_id = 'IFNG' THEN v.log2fc END) AS ifng_log2fc,
       MIN(p.status) AS status
FROM reference_values v JOIN reference_profiles p USING(profile_id)
GROUP BY p.cell_type_id, v.gene;
