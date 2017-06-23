### Mosaic SNV calling
This contains wrapper code to run Art's mosaic calling system. It does all the
messy bits of calling samtools, generating bcfs, fixing, merging and running
denovogear.

```sh
# grab the code
git clone https://github.com/jeremymcrae/mosaic_de_novos.git

# Then change directory and install the package:
cd mosaic_de_novos
python setup.py install

# Call mosaics for a single chromosome with:
python scripts/call_mosaics.py \
  --proband-bam PROBAND_BAM_PATH \
  --mother-bam MOTHER_BAM_PATH \
  --father-bam FATHER_BAM_PATH \
  --proband-sex male|female \
  --chrom 1 \
  --outdir RESULTS_DIR
```

You can also call variants in specific regions within a single chromosome with
`--start START_POS` and `--stop STOP_POS`. Otherwise it will default to the full
chromosome length. By default this looks for mosaic events with a proportion
around 0.25, but you can change this by including the argument `--proportion X`,
where X is between 0 and 1. NOTE: this might not be true, it might only use 0.25.

It should generate a file in RESULTS_DIR named:
`PROBAND_ID.denovogear.CHROM.START-STOP.modified.dnm`

Of course, this will probably fail in some way, most likely due to pysam
requirements. You can install pysam on the farm with:
```sh
pip install --user pysam
```

#### Merging
You can merge the denovogear output from different chromosomes with:
```sh
python scripts/merge_denovogear.py \
  --remove-files \ # removes the intermediate denovogear outputs
  --folder RESULTS_DIR \
  --pattern "modified"
```
