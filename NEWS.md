# tzara 0.0.9

* Added arguments `"seq_id"` and `"filename"` to `dadamap()`, to make it more
  user-friendly.  Use of `add_dada_names()` is now unnecessary.
* `extract_region()` auto-detects the format of the output file, if given,
  from the name.
* Fixed a bug where `extract_region()` was not saving quality information
  from the input to an output "fastq" file when the input was a "fastq" file or
  `QualityScaledXStringSet`.

# tzara 0.0.8

* Added a `NEWS.md` file to track changes to the package.
* **Breaking Change**: Uniformly use the column name "seq_id" rather than "seq",
  "seq.id", or "seq_name" for unique sequence identifiers.
* Add support for `DNAStringSet` and `QualityScaledDNAStringSet` in 
  `extract_region()`.
