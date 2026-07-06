# AlphaGenome ontology terms with RNA_SEQ tracks

Only ontology CURIEs that AlphaGenome actually has RNA_SEQ tracks for return
numeric values. Passing any other term → `Unable to extract numeric values for
'RNA_SEQ'`. This list is extracted from AlphaGenome's own example metadata
(output_type == RNA_SEQ). Regenerate/extend with:

    modal run modal_design.py::dump_ag_outputs

which prints `RNA_SEQ_ontology_to_biosample` (ontology_curie → biosample_name).

## Verified RNA_SEQ CL terms (ontology → biosample)
```
CL:0000047  neuronal stem cell
CL:0000062  osteoblast
CL:0000084  T-cell
CL:0000100  motor neuron
CL:0000115  endothelial cell
CL:0000121  Purkinje cell
CL:0000127  astrocyte
CL:0000134  mesenchymal stem cell
CL:0000137  osteocyte
CL:0000138  chondrocyte
CL:0000169  type B pancreatic cell
CL:0000182  hepatocyte                (HepG2)
CL:0000187  myocyte
CL:0000192  smooth muscle cell
CL:0000236  B cell
CL:0000312  keratinocyte
CL:0000623  natural killer cell
CL:0000624  CD4-positive, alpha-beta T cell
CL:0000625  CD8-positive, alpha-beta T cell
CL:0000746  cardiac muscle cell
CL:0000837  hematopoietic multipotent progenitor cell
CL:0000842  mononuclear cell
CL:0001053  IgD-negative memory B cell
CL:0001054  CD14-positive monocyte    (THP1)
CL:0001059  common myeloid progenitor, CD34-positive
CL:0002319  neural cell               (SH-SY5Y)
CL:0002327  mammary epithelial cell   (MCF7)
CL:0002518  kidney epithelial cell    (HEK293)
```
(Truncated; run `dump_ag_outputs` for the complete set — AlphaGenome also exposes
CAGE, ATAC, DNASE, CHIP_HISTONE, CHIP_TF, SPLICE_*, PROCAP, CONTACT_MAPS output types.)

## Cell-line → term mapping used in `elements.CELL_CONTEXTS`
| Cell line | ontology | note |
|-----------|----------|------|
| THP1  | CL:0001054 | CD14+ monocyte — good IFN/inflammation target |
| HepG2 | CL:0000182 | hepatocyte — glucocorticoid/oxidative/hypoxia |
| Jurkat| CL:0000084 | T-cell |
| MCF7  | CL:0002327 | mammary epithelial (ER+ context) |
| HEK293| CL:0002518 | kidney epithelial (generic) |
| SHSY5Y| CL:0002319 | neural cell |

Cell-line-specific EFO terms (e.g. HepG2=EFO:xxxx) may give sharper contrast if
present in the AlphaGenome build; confirm with `dump_ag_outputs` before using.
