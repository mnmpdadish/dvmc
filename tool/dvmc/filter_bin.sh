
# very basic script: to be upgraded

for f in output/zvo_nCHAm_nAHCm_*.bin; do
  convertOutputBin.py $f;
  dvmc_spectrum.py spectrumpara.def output 1;
done

