Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v0.0.1
---------------------------
+ Make specifying a region mandatory for the pipeline.
+ Use haplotype-shuffler instead of vcf-combo script.
+ Require bgzipped, tabix-indexed VCF files as input.
+ Add support for specifying a region to act on.
+ Add support for project configuration using
[PEP](http://pep.databio.org/en/latest/).
+ Add integration tests using pytest-workflow and github workflows.
