
### Mosaic SNV calling
This contains helper code to run Art's mosaic calling system. It does all the messy bits of calling samtools, generating bcfs, fixing, merging and running denovogear.

```sh
# grab the code
git clone https://github.com/jeremymcrae/mosaic_de_novos.git

# Then change directory:
cd mosaic_de_novos/src

# Call mosaics for a single chromosome with:
python mosaic_calling_denovo_gear.py \
  --proband-bam PROBAND_BAM_PATH \
  --mother-bam MOTHER_BAM_PATH \
  --father-bam FATHER_BAM_PATH \
  --proband-sex male \
  --chrom 1 \
  --outdir DIR_FOR_DENOVOGEAR_RESULTS
```

You can also call variants in specific regions within a single chromosome with `--start START_POS` and `--stop STOP_POS`. Otherwise it will default to the full chromosome length.

It should generate two files in DIR_FOR_DENOVOGEAR_RESULTS named:
`PROBAND_ID.denovogear.CHROM.START-STOP.standard.dnm`
`PROBAND_ID.denovogear.CHROM.START-STOP.modified.dnm`

Of course, this will probably fail in some way, most likely due to pysam requirements. You can install pysam on the farm with:
```sh
pip install --user pysam
```
