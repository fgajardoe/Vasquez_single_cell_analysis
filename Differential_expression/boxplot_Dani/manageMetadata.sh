zless GSE13204-GPL570_series_matrix.txt.gz |grep -ve'^!' > GSE13204-GPL570_series_matrix.noMetadata.tsv
zless GSE13204-GPL7473_series_matrix.txt.gz |grep -ve'^!' > GSE13204-GPL7473_series_matrix.noMetadata.tsv
zless GSE13204-GPL570_series_matrix.txt.gz |grep -e 'Sample_characteristics_ch1\|Sample_geo_accession' > GSE13204-GPL570_series_matrix.sampleMetadata.tsv
zless GSE13204-GPL7473_series_matrix.txt.gz |grep -e 'Sample_characteristics_ch1\|Sample_geo_accession' > GSE13204-GPL7473_series_matrix.sampleMetadata.tsv
