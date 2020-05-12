# tzara 0.0.9

* Added arguments `"seq_id"` and `"filename"` to `dadamap()`, to make it more
  user-friendly.  Use of `add_dada_names()` is now unnecessary.

# tzara 0.0.8

* Added a `NEWS.md` file to track changes to the package.
* **Breaking Change**: Uniformly use the column name "seq_id" rather than "seq",
  "seq.id", or "seq_name" for unique sequence identifiers.
* Add support for `DNAStringSet` and `QualityScaledDNAStringSet` in 
  `extract_region()`.