# Notes

- Include information on where this folder is, similar to what was done with `main_directory`.
- Set up prompt to commit series of versions, say of `MafA_SNP_DEG_Integration.R` from `v1` to `v6`.
  - Be careful with files that have implicit `v1`.

Find files ending in `[ _]v[0-9]+.R$` that have the same name before this version. If there is more than one set, do the following steps for each set.
Denote by `[basename]` the part of the filename before the version number.

1. If there is not a `v1` version but there is the `[basename]` without a version, copy `[basename].R` to `[basename] v1.R` or `[basename]_v1.R` as appropriate. If there is already a `[basename] v1.R` or `[basename]_v1.R`, do not overwrite it, but instead copy the `v1` version to `[basename].R` for this set.
2. Commit `[basename].R` with message `v1` and the date this file was created. (Be careful with files that have implicit `v1`.)
3. For subsequent versions, copy the next highest version number to `[basename].R` and commit with the next version number and the date this file was created. Continue until all versions are committed.
4. Do not include the singleton files such as `source*v2.R`.
5. Make sure to keep the versioned files intact.

- use `pandoc -s ex.doc -t markdown -o ex.md` to convert doc to markdown.
Find `docx` files and use pandoc to convert them to `md` versions
