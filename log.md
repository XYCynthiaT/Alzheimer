# April 26
- Add: glycolipid from Reactome
- Add: XXX_modified_reactome+glycolipids_(region).csv
- updata: regions.R, plots.R 

# April 30
- RNA.R: MSBB labels = :point_right: levels = 
- re-calculate everything related to MSBB

# May 9
- FIX: the rownames of logcpm didn't match rna-seq data
- - Update: xxx_normalized_adj.rds

# May 19 
- glyco_Reactome.R --> glyco.R, combing data from reactome and from review paper.
- Reactome: N-+O-+glycolipid = 462, 115 overlapped with glycosyltransferase list, 347 are not glycosyltransferases.
- Results in mayo_modified_paper_TCX.xlsx didn't match mayo_modified_reactome_glycolipids_TCX.xlsx --> update DE_mayo_adj.rds in app/data, update DE_rosmap_adj.rds in data and app/datas