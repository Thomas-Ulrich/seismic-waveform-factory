# dataset information

`dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest-fault.xdmf` is obtained
 through the `rapid-earthquake-dynamics`. The folder name (with earthquake
 information) is `2024-01-22_Mw7_China_us7000lsze_1`.
The data is extracted with :

```python
seissol_output_extractor \
extracted_output/dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extracted-fault.xdmf \
--var ASl Sls Sld  --time "i::8" --add test
```

Note that the dataset is subsampled in time to reduce its size, and is meant
 only for testing for the workflow.
